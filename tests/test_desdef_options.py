from msd_interface_design import DesignDefinitionOpts
from opt_holder import OptHolder
import blargs


def test_create_desdef_options():
    d = {}
    p = blargs.Parser(d)
    DesignDefinitionOpts.add_options(p)
    assert "des_def" in p._options
    p.process_command_line(["--des_def", "foo"])
    assert "des_def" in d
    assert d["des_def"] == "foo"


def test_optholder():
    opts = OptHolder()
    p = blargs.Parser(opts)
    p.int("foo")
    p.process_command_line(["--foo", "12"])
    assert hasattr(opts, "foo")
    assert opts.foo == 12
