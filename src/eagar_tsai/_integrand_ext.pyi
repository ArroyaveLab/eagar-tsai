"""Type stub for the _integrand_ext C extension module."""

def get_integrand_capsule() -> object:
    """Return a PyCapsule wrapping the Eagar-Tsai integrand.

    The capsule's name is ``"double (int, double *)"`` — one of the four
    signatures recognised by :class:`scipy.integrate.LowLevelCallable`.

    Returns:
        A PyCapsule suitable for wrapping with ``LowLevelCallable``.
    """
    ...
