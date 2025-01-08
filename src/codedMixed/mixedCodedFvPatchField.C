/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mixedCodedFvPatchField.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void mixedCodedFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", redirectType_);
    setFieldTemplates(dynCode);
    //dynCode.setFilterVariable("codeWrite", codeWrite_);

    // compile filtered C template
    dynCode.addCompileFile("mixedCodedBCTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("mixedCodedBCTemplate.H");
    //dynCode.addCopyFile("mixedCodedFvPatchFieldsFwdTemplate.H");
    //dynCode.addCopyFile("mixedCodedFvPatchFieldsTemplate.H");
    //dynCode.addCopyFile("mixedCodedFvPatchFieldTemplate.H");
    //dynCode.addCopyFile("mixedCodedFvPatchFieldTemplate.C");

    // debugging: make BC verbose
    //         dynCode.setFilterVariable("verbose", "true");
    //         Info<<"compile " << redirectType_ << " sha1: "
    //             << context.sha1() << endl;

    // define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
            "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lfoam \\\n"
            + "    -lfiniteVolume \\\n"
            + "    -lmeshTools \\\n"
            + context.libs()
        );
}


template <class Type>
Foam::dlLibraryTable& mixedCodedFvPatchField<Type>::libs() const
{
    return const_cast<Time&>(this->db().time()).libs();
}


template <class Type>
Foam::string mixedCodedFvPatchField<Type>::description() const
{
    return "mixedCodeBC referring to " + redirectType_;
}


template <class Type>
void mixedCodedFvPatchField<Type>::clearRedirect() const
{
    redirectBCPtr_.clear();
}


template <class Type>
const Foam::dictionary& mixedCodedFvPatchField<Type>::codeDict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
mixedCodedFvPatchField<Type>::mixedCodedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    codedBase(),
    dict_(),
    redirectType_()
{}


template<class Type>
mixedCodedFvPatchField<Type>::mixedCodedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    codedBase(),
    dict_(dict),
    redirectType_(dict.lookup("redirectType"))
{
    // Call evaluate only if the value is not found. Used to avoid evaluating
    // when we have incomplete meshes during Parallel Load Balancing. When
    // shipping the field over to another processor, we first call write, making
    // sure that the value is written and read it on the other side (see
    // write member function). If this proves to be problematic, we can always
    // initialize with patch internal field for the start-up. VV, 12/Apr/2019.
    if (!dict.found("value"))
    {
        evaluate();
    }
}


template<class Type>
mixedCodedFvPatchField<Type>::mixedCodedFvPatchField
(
    const mixedCodedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_)
{}


template<class Type>
mixedCodedFvPatchField<Type>::mixedCodedFvPatchField
(
    const mixedCodedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_)
{}


template<class Type>
mixedCodedFvPatchField<Type>::mixedCodedFvPatchField
(
    const mixedCodedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixedCodedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void mixedCodedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
void mixedCodedFvPatchField<Type>::evaluate(const Pstream::commsTypes comms)
{
    updateLibrary(redirectType_);
    return redirectBC().evaluate(comms);
}


template<class Type>
tmp<Field<Type> > mixedCodedFvPatchField<Type>::snGrad() const
{
    updateLibrary(redirectType_);
    return redirectBC().snGrad();
}


template<class Type>
tmp<Field<Type> > mixedCodedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& tf
) const
{
    updateLibrary(redirectType_);
    return redirectBC().valueInternalCoeffs(tf);
}


template<class Type>
tmp<Field<Type> > mixedCodedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& tf
) const
{
    updateLibrary(redirectType_);
    return redirectBC().valueBoundaryCoeffs(tf);
}


template<class Type>
tmp<Field<Type> > mixedCodedFvPatchField<Type>::gradientInternalCoeffs() const
{
    updateLibrary(redirectType_);
    return redirectBC().gradientInternalCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedCodedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    updateLibrary(redirectType_);
    return redirectBC().gradientBoundaryCoeffs();
}


template<class Type>
void mixedCodedFvPatchField<Type>::write(Ostream& os) const
{
    //mixedFvPatchField::write(os);
    os << dict_ << endl;
}

template<class Type>
fvPatchField<Type>& mixedCodedFvPatchField<Type>::redirectBC() const
{
    if (!redirectBCPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", redirectType_);

        redirectBCPtr_ = fvPatchField<Type>::New
        (
            redirectType_,
            this->patch(),
            this->dimensionedInternalField()
        );
    }
    return redirectBCPtr_();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
