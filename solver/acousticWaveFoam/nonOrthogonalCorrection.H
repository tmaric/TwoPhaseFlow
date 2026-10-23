// Explicit non-orthogonal correction of the pressure Laplacian and its
// transpose, for the monolithic update.
//
// A corrected Laplacian splits the volume-integrated operator into the
// orthogonal part K, assembled in the block matrix, and the correction
//   C p = D B I G p,
// with G the Gauss cell gradient, I the linear interpolation to the faces,
// B the face flux gamma |Sf| k_f . g_f with the non-orthogonal correction
// vectors k_f, and D the sum of face fluxes over each cell. C equals the
// correction of the Gauss linear corrected scheme. C is not symmetric, so
// K + C has complex eigenvalues on non-orthogonal meshes, and the lossless
// wave equation then has growing modes. (C + C^T)/2 keeps K + (C + C^T)/2
// symmetric. The transpose is formed from the transposes of the four
// factors, with the linear part of the boundary values (the boundary value
// follows the cell value except on fixed-value patches).
//
// Processor patches are faces between two cells: each processor adds the
// contribution of its own cell, with the neighbour values (cell values,
// gradients, scaled gradients) exchanged through temporary fields. Other
// coupled patches (cyclic, AMI) are not supported.

#include "emptyFvPatch.H"
#include "processorFvPatch.H"
#include "calculatedFvPatchFields.H"

// Cell values with boundary values of p's patch types: the actual boundary
// values of p, or for the linear part zero on fixed-value patches.
static tmp<volScalarField> nonOrthScalarField
(
    const volScalarField& p,
    const scalarField& values,
    const bool linearPart
)
{
    const fvMesh& mesh = p.mesh();
    const IOobject io
    (
        "nonOrthScalar",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    tmp<volScalarField> tf;
    if (linearPart)
    {
        wordList types(p.boundaryField().types());
        forAll(types, patchi)
        {
            if (p.boundaryField()[patchi].fixesValue())
            {
                types[patchi] = "fixedValue";
            }
        }
        tf = tmp<volScalarField>
        (
            new volScalarField
            (
                io, mesh, dimensionedScalar(p.dimensions(), Zero), types
            )
        );
    }
    else
    {
        tf = tmp<volScalarField>(new volScalarField(io, p));
    }
    tf.ref().primitiveFieldRef() = values;
    tf.ref().correctBoundaryConditions();
    return tf;
}


// Cell vectors with their neighbour values on the processor patches.
static tmp<volVectorField> nonOrthVectorHalo
(
    const fvMesh& mesh,
    const vectorField& values
)
{
    tmp<volVectorField> tf
    (
        new volVectorField
        (
            IOobject
            (
                "nonOrthVector",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector(dimless, Zero),
            calculatedFvPatchVectorField::typeName
        )
    );
    tf.ref().primitiveFieldRef() = values;
    tf.ref().correctBoundaryConditions();
    return tf;
}


static void nonOrthCheckPatches(const fvMesh& mesh)
{
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (patch.coupled() && !isA<processorFvPatch>(patch))
        {
            FatalErrorInFunction
                << "nonOrthogonalCorrection supports processor patches only, "
                << "not the coupled patch " << patch.name()
                << exit(FatalError);
        }
    }
}


// Gauss gradient of a field whose boundary values have been set.
static vectorField nonOrthGaussGradient(const volScalarField& f)
{
    const fvMesh& mesh = f.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();
    const scalarField& values = f.primitiveField();

    vectorField sum(mesh.nCells(), Zero);

    forAll(owner, facei)
    {
        const scalar w = weights[facei];
        const vector flux =
            Sf[facei]
           *(w*values[owner[facei]] + (1 - w)*values[neighbour[facei]]);
        sum[owner[facei]] += flux;
        sum[neighbour[facei]] -= flux;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch))
        {
            continue;
        }
        const fvPatchScalarField& pf = f.boundaryField()[patchi];
        const labelUList& cells = patch.faceCells();
        const vectorField& Sb = Sf.boundaryField()[patchi];

        if (patch.coupled())
        {
            const scalarField& w = weights.boundaryField()[patchi];
            const tmp<scalarField> tnbr = pf.patchNeighbourField();
            const scalarField& nbr = tnbr();
            forAll(cells, facei)
            {
                sum[cells[facei]] +=
                    Sb[facei]
                   *(w[facei]*values[cells[facei]] + (1 - w[facei])*nbr[facei]);
            }
        }
        else
        {
            forAll(cells, facei)
            {
                sum[cells[facei]] += Sb[facei]*pf[facei];
            }
        }
    }

    return sum/mesh.V().field();
}


// Face fluxes gamma |Sf| k_f . g_f of the correction, on the internal faces
// and on the processor patches (local orientation); zero elsewhere.
static tmp<surfaceScalarField> nonOrthCorrectionFlux
(
    const volScalarField& f,
    const surfaceScalarField& gammaMagSf
)
{
    const fvMesh& mesh = f.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& corrVecs = mesh.nonOrthCorrectionVectors();

    const tmp<volVectorField> tgrad
    (
        nonOrthVectorHalo(mesh, nonOrthGaussGradient(f))
    );
    const volVectorField& grad = tgrad();

    tmp<surfaceScalarField> tflux
    (
        new surfaceScalarField
        (
            IOobject
            (
                "nonOrthCorrectionFlux",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar
            (
                gammaMagSf.dimensions()*f.dimensions()/dimLength,
                Zero
            )
        )
    );
    surfaceScalarField& flux = tflux.ref();

    forAll(owner, facei)
    {
        const scalar w = weights[facei];
        flux[facei] =
            gammaMagSf[facei]
           *(corrVecs[facei]
           & (w*grad[owner[facei]] + (1 - w)*grad[neighbour[facei]]));
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (!patch.coupled())
        {
            continue;
        }
        const scalarField& w = weights.boundaryField()[patchi];
        const vectorField& k = corrVecs.boundaryField()[patchi];
        const scalarField& gm = gammaMagSf.boundaryField()[patchi];
        const tmp<vectorField> tnbr =
            grad.boundaryField()[patchi].patchNeighbourField();
        const vectorField& nbr = tnbr();
        const labelUList& cells = patch.faceCells();
        scalarField& pflux = flux.boundaryFieldRef()[patchi];

        forAll(cells, facei)
        {
            pflux[facei] =
                gm[facei]
               *(k[facei]
               & (w[facei]*grad[cells[facei]] + (1 - w[facei])*nbr[facei]));
        }
    }

    return tflux;
}


// Sum of the face fluxes over each cell (D).
static scalarField nonOrthFaceSum(const surfaceScalarField& flux)
{
    const fvMesh& mesh = flux.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    scalarField result(mesh.nCells(), Zero);
    forAll(owner, facei)
    {
        result[owner[facei]] += flux[facei];
        result[neighbour[facei]] -= flux[facei];
    }
    forAll(mesh.boundary(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            continue;
        }
        const labelUList& cells = mesh.boundary()[patchi].faceCells();
        const scalarField& pflux = flux.boundaryField()[patchi];
        forAll(cells, facei)
        {
            result[cells[facei]] += pflux[facei];
        }
    }
    return result;
}


// C^T phi for the linear part of C.
static scalarField nonOrthTranspose
(
    const volScalarField& p,
    const scalarField& values,
    const surfaceScalarField& gammaMagSf
)
{
    const fvMesh& mesh = p.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& corrVecs = mesh.nonOrthCorrectionVectors();
    const surfaceVectorField& Sf = mesh.Sf();

    // D^T and B^T give a face vector h; I^T distributes it to the cells.
    const tmp<volScalarField> tphi(nonOrthScalarField(p, values, true));
    const volScalarField& phi = tphi();

    vectorField X(mesh.nCells(), Zero);
    forAll(owner, facei)
    {
        const vector h =
            gammaMagSf[facei]*corrVecs[facei]
           *(values[owner[facei]] - values[neighbour[facei]]);
        const scalar w = weights[facei];
        X[owner[facei]] += w*h;
        X[neighbour[facei]] += (1 - w)*h;
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (!patch.coupled())
        {
            continue;
        }
        const scalarField& w = weights.boundaryField()[patchi];
        const vectorField& k = corrVecs.boundaryField()[patchi];
        const scalarField& gm = gammaMagSf.boundaryField()[patchi];
        const tmp<scalarField> tnbr =
            phi.boundaryField()[patchi].patchNeighbourField();
        const scalarField& nbr = tnbr();
        const labelUList& cells = patch.faceCells();
        forAll(cells, facei)
        {
            X[cells[facei]] +=
                w[facei]*gm[facei]*k[facei]*(values[cells[facei]] - nbr[facei]);
        }
    }

    // G^T = Gs^T V^{-1}.
    const tmp<volVectorField> tY
    (
        nonOrthVectorHalo(mesh, X/mesh.V().field())
    );
    const volVectorField& Y = tY();

    scalarField result(mesh.nCells(), Zero);
    forAll(owner, facei)
    {
        const scalar w = weights[facei];
        const scalar s =
            (Y[owner[facei]] - Y[neighbour[facei]]) & Sf[facei];
        result[owner[facei]] += w*s;
        result[neighbour[facei]] += (1 - w)*s;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if (isA<emptyFvPatch>(patch))
        {
            continue;
        }
        const labelUList& cells = patch.faceCells();
        const vectorField& Sb = Sf.boundaryField()[patchi];

        if (patch.coupled())
        {
            const scalarField& w = weights.boundaryField()[patchi];
            const tmp<vectorField> tnbr =
                Y.boundaryField()[patchi].patchNeighbourField();
            const vectorField& nbr = tnbr();
            forAll(cells, facei)
            {
                result[cells[facei]] +=
                    w[facei]*((Y[cells[facei]] - nbr[facei]) & Sb[facei]);
            }
        }
        else if (!p.boundaryField()[patchi].fixesValue())
        {
            forAll(cells, facei)
            {
                result[cells[facei]] += Y[cells[facei]] & Sb[facei];
            }
        }
    }

    return result;
}


// The correction used in the q equation: C p (explicit) or
// (C p + C^T p)/2 (symmetric), where the transpose acts on the linear part.
// p must have up-to-date boundary values.
static scalarField nonOrthCorrection
(
    const volScalarField& p,
    const surfaceScalarField& gammaMagSf,
    const word& mode
)
{
    const scalarField Cp
    (
        nonOrthFaceSum(nonOrthCorrectionFlux(p, gammaMagSf)())
    );
    if (mode == "explicit")
    {
        return Cp;
    }
    return 0.5*(Cp + nonOrthTranspose(p, p.primitiveField(), gammaMagSf));
}


// Adjoint test <phi, C psi> = <C^T phi, psi> for random fields and the
// linear part of C; returns the relative difference.
static scalar nonOrthAdjointTest
(
    const volScalarField& p,
    const surfaceScalarField& gammaMagSf
)
{
    const fvMesh& mesh = p.mesh();
    Random rnd(1234 + Pstream::myProcNo());
    scalarField psi(mesh.nCells()), phi(mesh.nCells());
    forAll(psi, celli)
    {
        psi[celli] = rnd.sample01<scalar>() - 0.5;
        phi[celli] = rnd.sample01<scalar>() - 0.5;
    }

    const scalarField Cpsi
    (
        nonOrthFaceSum
        (
            nonOrthCorrectionFlux
            (
                nonOrthScalarField(p, psi, true)(),
                gammaMagSf
            )()
        )
    );
    const scalarField CTphi(nonOrthTranspose(p, phi, gammaMagSf));
    const scalar lhs = gSum(phi*Cpsi);
    const scalar rhs = gSum(CTphi*psi);
    return mag(lhs - rhs)/max(max(mag(lhs), mag(rhs)), VSMALL);
}
