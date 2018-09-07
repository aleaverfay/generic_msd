class OptHolder:
    """A simple class to use with blargs

    Turns every option that the blargs.Parser class defines into
    an attribute of the intsance. It does this by interpretting
    __setitem__ as setattr
    """

    def __setitem__(self, key, val):
        setattr(self, key, val)
