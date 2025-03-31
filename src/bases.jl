BaseTypes = (
    linearLagrange = "Lag1",
    quadraticLagrange = "Lag2",
    cubicLagrange = "Lag3",
    cubicHermite = "Her3",
)

LocalBases = (
    Lag1 = (ne) -> (type = BaseTypes.linearLagrange, p = 1, nB = 2, neq = ne - 1),
    Lag2 = (ne) -> (type = BaseTypes.quadraticLagrange, p = 2, nB = 3, neq = 2*ne - 1),
    Lag3 = (ne) -> (type = BaseTypes.cubicLagrange, p = 3, nB = 4, neq = 3*ne - 1),
    Her3 = (ne) -> (type = BaseTypes.cubicHermite, p = 3, nB = 4, neq = 2*ne)
)

function monta_base(baseType, ne)
    return LocalBases[Symbol(baseType)](ne)
end