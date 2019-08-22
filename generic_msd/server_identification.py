# What computer am I running on?!

from enum import Enum
import socket
import toolz.functoolz
import traceback
import sys


class KnownComputers(Enum):
    KILLDEVIL = "killdevil"
    LONGLEAF = "longleaf"
    DOGWOOD = "dogwood"
    WIGGINS = "wiggins"
    ANDREWS_LAPTOP = "andrews_laptop"
    CATHYS_DESKTOP = "cathys_desktop"

class ServerIdentifier:
    def __init__(self, *, masq=None):
        self.masquerade = masq

    @toolz.functoolz.memoize
    def what_computer(self) -> KnownComputers:
        """Return the id for the computer this script is running on.

        This function is ever so slightly expensive, and needs only to
        be run once, so it is memoized"""
        if self.masquerade:
            return self.masquerade

        hostname = socket.gethostname()
        if self._on_killdevil(hostname):
            print("on killdevil")
            return KnownComputers.KILLDEVIL
        elif self._on_dogwood(hostname):
            print("on dogwood")
            return KnownComputers.DOGWOOD
        elif hostname == "wiggins":
            print("on wiggins")
            return KnownComputers.WIGGINS
        elif hostname.startswith("Lysis"):
            print("on lysis")
            # we're on andrew's laptop
            return KnownComputers.ANDREWS_LAPTOP
        else:
            # assume we're on Cathy's desktop?
            print("what computer? Assume Cathy's desktop", hostname)
            return KnownComputers.CATHYS_DESKTOP

    def _on_killdevil(self, hostname=None):
        """Logic for figuring out whether or not you're on killdevil vs dogwood"""
        if hostname is None:
            hostname = socket.gethostname()[:9]
        if hostname.startswith("killdevil"):
            return True
        elif hostname.startswith("c-"):
            # ok, distinguish between dogwood and killdevil:
            node_num = int(hostname.split("-")[1])
            return node_num >= 183 and node_num <= 199
        else:
            return False

    def _on_dogwood(self, hostname=None):
        """Logic for figuring out whether or not you're on dogwood vs killdevil"""
        if hostname is None:
            hostname = socket.gethostname()
        if hostname.startswith("dogwood"):
            return True
        elif hostname.startswith("c-"):
            node_num = int(hostname.split("-")[1])
            return node_num >= 201 and node_num <= 211
        else:
            return False
