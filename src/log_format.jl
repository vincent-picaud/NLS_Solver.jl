#
# Define how to format logger outputs
#

using Printf: @sprintf

_fmt(x::Integer) = @sprintf("%4d", x)
_fmt(x::AbstractFloat) = @sprintf("% .6e", x)
