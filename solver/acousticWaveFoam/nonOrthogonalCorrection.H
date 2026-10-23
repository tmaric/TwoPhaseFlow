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
// Serial meshes without coupled patches only.

#include "emptyFvPatch.H"

static vectorField nonOrthGaussGradient
(
    const volScalarField& p,
    const scalarField& values,
    const bool linearPart
)
{
    const fvMesh& mesh = p.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();

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
        const fvPatchScalarField& pp = p.boundaryField()[patchi];
        const labelUList& cells = patch.faceCells();
        const vectorField& Sb = Sf.boundaryField()[patchi];

        forAll(cells, facei)
        {
            const scalar value =
                linearPart
              ? (pp.fixesValue() ? 0 : values[cells[facei]])
              : pp[facei];
            sum[cells[facei]] += Sb[facei]*value;
        }
    }

    return sum/mesh.V().field();
}


// Face fluxes gamma |Sf| k_f . g_f of the correction on the internal faces.
static scalarField nonOrthCorrectionFlux
(
    const volScalarField& p,
    const surfaceScalarField& gammaMagSf
)
{
    const fvMesh& mesh = p.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& corrVecs = mesh.nonOrthCorrectionVectors();
    const vectorField grad(nonOrthGaussGradient(p, p.primitiveField(), false));

    scalarField flux(owner.size());
    forAll(owner, facei)
    {
        const scalar w = weights[facei];
        flux[facei] =
            gammaMagSf[facei]
           *(corrVecs[facei]
           & (w*grad[owner[facei]] + (1 - w)*grad[neighbour[facei]]));
    }
    return flux;
}


// C p (transpose = false) or C^T phi (transpose = true).
static scalarField applyNonOrthCorrection
(
    const volScalarField& p,
    const scalarField& values,
    const surfaceScalarField& gammaMagSf,
    const bool transpose
)
{
    const fvMesh& mesh = p.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& corrVecs = mesh.nonOrthCorrectionVectors();
    const surfaceVectorField& Sf = mesh.Sf();

    scalarField result(mesh.nCells(), Zero);

    if (!transpose)
    {
        const vectorField grad(nonOrthGaussGradient(p, values, false));
        forAll(owner, facei)
        {
            const scalar w = weights[facei];
            const vector g =
                w*grad[owner[facei]] + (1 - w)*grad[neighbour[facei]];
            const scalar flux =
                gammaMagSf[facei]*(corrVecs[facei] & g);
            result[owner[facei]] += flux;
            result[neighbour[facei]] -= flux;
        }
        return result;
    }

    // Transpose: D^T, B^T, I^T, then G^T = Gs^T V^{-1}.
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
    const vectorField Y(X/mesh.V().field());

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
        if (isA<emptyFvPatch>(patch) || p.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& cells = patch.faceCells();
        const vectorField& Sb = Sf.boundaryField()[patchi];
        forAll(cells, facei)
        {
            result[cells[facei]] += Y[cells[facei]] & Sb[facei];
        }
    }

    return result;
}


// The correction used in the q equation: C p (explicit) or
// (C p + C^T p)/2 (symmetric), where the transpose acts on the linear part.
static scalarField nonOrthCorrection
(
    const volScalarField& p,
    const surfaceScalarField& gammaMagSf,
    const word& mode
)
{
    const scalarField Cp
    (
        applyNonOrthCorrection(p, p.primitiveField(), gammaMagSf, false)
    );
    if (mode == "explicit")
    {
        return Cp;
    }
    return 0.5*
    (
        Cp
      + applyNonOrthCorrection(p, p.primitiveField(), gammaMagSf, true)
    );
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
    Random rnd(1234);
    scalarField psi(mesh.nCells()), phi(mesh.nCells());
    forAll(psi, celli)
    {
        psi[celli] = rnd.sample01<scalar>() - 0.5;
        phi[celli] = rnd.sample01<scalar>() - 0.5;
    }

    // Linear part of C psi: boundary values follow the cells.
    const vectorField grad(nonOrthGaussGradient(p, psi, true));
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();
    const surfaceVectorField& corrVecs = mesh.nonOrthCorrectionVectors();
    scalarField Cpsi(mesh.nCells(), Zero);
    forAll(owner, facei)
    {
        const scalar w = weights[facei];
        const scalar flux =
            gammaMagSf[facei]
           *(corrVecs[facei]
           & (w*grad[owner[facei]] + (1 - w)*grad[neighbour[facei]]));
        Cpsi[owner[facei]] += flux;
        Cpsi[neighbour[facei]] -= flux;
    }

    const scalarField CTphi
    (
        applyNonOrthCorrection(p, phi, gammaMagSf, true)
    );
    const scalar lhs = sum(phi*Cpsi);
    const scalar rhs = sum(CTphi*psi);
    return mag(lhs - rhs)/max(max(mag(lhs), mag(rhs)), VSMALL);
}
