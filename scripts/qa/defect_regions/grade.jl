# Grading rules applied to the raw sweep output.  Shared by analyze.jl and
# make_plots.jl so the report and the figures can never disagree.

const CONF = 6.0        # a product "confirms" a region at |S| >= CONF
const MARG = 3.0        # "marginal" support
const MAXEXT = 256      # a defect *region* is compact: bbox <= MAXEXT on each axis
const MINFILL = 0.25    # ... and fills at least this fraction of its bounding box

const PRODUCTS = ("dark_rate", "flat_im", "flat2d")

ev(r, nm) = (i = findfirst(e -> e.name == nm, r.evidence); i === nothing ? nothing : r.evidence[i])
sig(r, nm) = (e = ev(r, nm); e === nothing ? 0.0 : e.maxS)
medz(r, nm) = (e = ev(r, nm); e === nothing ? 0.0 : e.medz)
nconf(r) = count(nm -> abs(sig(r, nm)) >= CONF, PRODUCTS)
nmarg(r) = count(nm -> MARG <= abs(sig(r, nm)) < CONF, PRODUCTS)

"""Corroboration tier.  A+ : all three products confirm; A : two of three;
B : one confirms and at least one other is marginal; C : a single product only."""
tier(r) = nconf(r) >= 3 ? "A+" : nconf(r) >= 2 ? "A" :
          (nconf(r) == 1 && nmarg(r) >= 1) ? "B" : "C"

wx(r) = r.x1 - r.x0 + 1
wy(r) = r.y1 - r.y0 + 1
fillfrac(r) = r.area / (wx(r) * wy(r))

"""A defect *region* is compact.  Extended coherent structure — edge roll-off, broad
dark-current gradients — is real and the finder does detect it, but it is what the flat
field and dark subtraction exist to handle, not something to put in a bad-pixel mask;
it is reported separately rather than silently dropped."""
compact(r) = wx(r) <= MAXEXT && wy(r) <= MAXEXT && fillfrac(r) >= MINFILL

"""The set proposed to AKS as a mask."""
proposed(r) = tier(r) in ("A", "A+") && compact(r)

severity(r) = r.area * maximum(abs(sig(r, n)) for n in PRODUCTS)
