"""Foreground/background HI columns and optical depths from the Soeding et al.
(2025) 3D reconstruction of the Milky Way.

Basic use::

    import hi3dtau

    hi3dtau.p(120.0, 30.0, 1.0)
    hi3dtau.columns(120.0, 30.0, 1.0)
    hi3dtau.tau(120.0, 30.0, 1.0, velocity_kms=0.0)

Data products: Soeding, Edenhofer, Ensslin et al. (2025), A&A 693, A139,
Zenodo doi:10.5281/zenodo.12578443.
"""

from .api import HI3D, columns, get, p, summary, tau
from .data import download

__all__ = ["HI3D", "columns", "download", "get", "p", "summary", "tau"]
__version__ = "0.1.0"
