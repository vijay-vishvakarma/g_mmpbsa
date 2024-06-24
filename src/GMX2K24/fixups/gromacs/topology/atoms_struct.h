typedef struct t_atom
{
    real           m, q;       /* Mass and charge                      */
    real           mB, qB;     /* Mass and charge for Free Energy calc */
    unsigned short type;       /* Atom type                            */
    unsigned short typeB;      /* Atom type for Free Energy calc       */
    ParticleType   ptype;      /* Particle type                        */
    int            resind;     /* Index into resinfo (in t_atoms)      */
    int            atomnumber; /* Atomic Number or 0                   */
    char           elem[4];    /* Element name                         */
    int            nrexcl;     /* Number of exclusions for this atom   */
    int*           excl;       /* List of excluded atom indices        */
} t_atom;