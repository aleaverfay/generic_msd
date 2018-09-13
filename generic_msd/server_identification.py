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
    @toolz.functoolz.memoize
    def what_computer(self) -> KnownComputers:
        """Return the id for the computer this script is running on.

        This function is ever so slightly expensive, and needs only to
        be run once, so it is memoized"""
        hostname = socket.gethostname()
        if self._on_killdevil(hostname):
            return KnownComputers.KILLDEVIL
        elif self._on_dogwood(hostname):
            return KnownComputers.DOGWOOD
        elif hostname == "wiggins":
            return KnownComputers.WIGGINS
        elif hostname.startswith("Lysis"):
            # we're on andrew's laptop
            return KnownComputers.ANDREWS_LAPTOP
        else:
            # assume we're on Cathy's desktop?
            print("what computer? Assume Cathy's desktop", hostname)
            return KnownComputers.CATHYS_DESKTOP

    def _on_killdevil(hostname=None):
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

    def _on_dogwood(hostname=None):
        """Logic for figuring out whether or not you're on dogwood vs killdevil"""
        if hostname is None:
            hostname = socket.gethostname()
        if hostname.startswith("dogwood"):
            return True
        elif hostname.startswith("c-"):
            node_num = int(hostname.split("-")[1])
            return node_num >= 201 and node_num <= 209
        else:
            return False
