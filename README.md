# [WIP] Coded Boundary Condition for Foam-Extend 5.0

This is simply a coded BC based on mixed BC in Foam-extend 5.0

> [!Note]
> Not well tested, PoC at best - Needs more actual testing

## Setup and usage

### Compile

```bash
source <your-foam-extend-5-bashrc>
wmake libso src/codedMixed
```

This compiles a small library `libcodedBC.so` into your `$FOAM_USER_LIBBIN`

## Run a case with the custom coded BC

In your shell:
```bash
export FOAM_CODE_TEMPLATES=<path-to-this-repo>/etc/codeTemplates/mixedCodedTemplate
```

In your case's `system/controlDict`
```cpp
libs ( "libcodedBC.so" );
InfoSwitches
{
    allowSystemOperations 1;
}
```

In the field file, you use the BC as follows:
```cpp
patchName
{
    type            mixedCoded;   // coded BC type name
    redirectType    custom;       // arbitrary type name so you can have multiples
    code #{
        // If this is for a scalar field, the patch data are scalars
        // Same if it's a vector or tensor field, the patch data automatically follow
        // the data type of the field this BC is attached to
        Info << " --- ---- ---- Executing custom BC" << endl;
        // But gradients are always vectors I guess
        // This mimics "zeroGradient"
        refGrad() = Field<vector>(this->patch().size(), vector::zero);
    #};  
}
```

There is a [tutorials/cavity](tutorials/cavity) case to showcase how this works:
```bash
# from tutorials/cavity
export FOAM_CODE_TEMPLATES=$PWD/../../etc/codeTemplates/mixedCodedTemplate
blockMesh && icoFoam
```
