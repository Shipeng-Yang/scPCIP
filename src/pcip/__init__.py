"""scPCIP: program-based cell-type annotation for plant single-cell RNA-seq.

Cell identity is scored with weighted gene programs over cross-species orthogroups rather than single marker genes.
"""

__version__ = "1.0.0"
__all__ = ["__version__", "annotate"]


def _greet() -> None:
    import sys

    from . import banner
    if banner.enabled() and banner.interactive_session():
        banner.show(sys.stdout)


_greet()


def annotate(*args, **kwargs):
    """See :func:`pcip.api.annotate`."""
    from .api import annotate as _annotate

    return _annotate(*args, **kwargs)

