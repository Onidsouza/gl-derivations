"""

Exact rational computations with gl(n) and its symmetric powers.

Provides Lie algebra elements, symmetric powers with the induced adjoint
action, and subspaces for computing invariants and weight decompositions.

"""

from .lie import GeneralLinear, LieElement
from .symmetric import SymmetricPower, SymmetricElement
from .subspaces import Subspace
from .actions import act, action_matrix, evaluate

__all__ = [
    "GeneralLinear",
    "LieElement",
    "SymmetricPower",
    "SymmetricElement",
    "Subspace",
    "act",
    "action_matrix",
    "evaluate"
]