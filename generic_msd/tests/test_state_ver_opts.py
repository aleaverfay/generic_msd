from msd_interface_design import StateVersionOpts
from opt_holder import OptHolder
import blargs


def test_state_ver_opts():
    opts = OptHolder()
    p = blargs.Parser(opts)
    StateVersionOpts.add_options(p)
    assert "state_version" in p._options
    p.process_command_line(["--state_version", "bar"])
    assert hasattr(opts, "state_version")
    assert opts.state_version == "bar"

    p.process_command_line(["-s", "baz"])
    assert opts.state_version == "baz"
