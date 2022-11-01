"""
Module for generating, writing, and reading opaque blinding factors
"""
from array import array
import numpy as np
import argparse
import os
import sha256

np.seterr(over='ignore')

def main():
    """
    Generates a blinding factor, saves it to disk, and verifies that reading was successful.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("seed_str", type=str, help="the seed string")
    seed_str = parser.parse_args().seed_str

    fname = f"{seed_str}.bin"
    blind = compute_blinding_factor(seed_str)
    write_blind(fname, blind)
    blind2 = read_blind(fname)

    print("Created blind =", blind)
    print("Wrote blind to", fname)
    print("Read blind from", fname)
    print("Read blind =", blind2)
    print("Blinds matched before and after reading?", blind == blind2)
    print("Blind between 0.95 and 1.05?", (0.95 < blind) and (blind < 1.05))


class Blind(float):
    """
    An "opaque" version of float which hides its value in print statements but can be used normally
    in arithematic.
    """
    def __new__(cls, value):
        return super().__new__(cls, value)

    def __init__(self, value):
        float.__init__(value)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "Blind(<opaque value>)"


def seed_from_string(seed_str):
    """
    Creates a deterministically random seed from an input string using hashing.
    Credit: StackExchange user jakevdp
    https://stackoverflow.com/questions/36755940/seed-numpy-random-randomstate-with-hashlib-hash
    Args:
        seed_str: str, the string to hash and convert to the seed
    Returns:
        seed: array of integers to be used with np.random.seed(seed)
    """
    assert isinstance(seed_str, str)
    ahash = sha256.sha256(seed_str)
    seed = np.frombuffer(ahash.hexdigest().encode(), dtype='uint32')[0]
    print(seed)
    return seed


def compute_blinding_factor(seed_str):
    """
    Computes a deterministically random blinding factor between 0.95 and 1.05 given an input string,
    which is hashed to seed the random number generator.
    Args:
        seed_str: str, the string to hash and use as the seed
    Returns:
        blind: Blind between 0.95 and 1.05
    """
    seed = seed_from_string(seed_str)
    png = PRN(seed, 1)
    return Blind(png.uniform(0.95, 1.05))


def write_blind(fname, blind):
    """
    Writes the blinding factor (in binary) to file.
    Args:
        fname: str, the full path to the file
        blind: float, the blinding factor
    """
    if os.path.isfile(fname):
        raise FileExistsError(f"{fname} already exists.")
    with open(fname, 'wb') as ofile:
        float_array = array('d', [blind])
        float_array.tofile(ofile)


def read_blind(fname):
    """
    Reads the blinding factor from file.
    Args:
        fname: str, the full path to the file
    Returns:
        Blind
    """
    with open(fname, 'rb') as ifile:
        blind = array('d')
        blind.frombytes(ifile.read())
        blind, = blind
    return Blind(blind)

class PRN:
    """
    Pseudorandom number generator.
    This code is ported, with only minor revision from the MILC code,
    https://github.com/milc-qcd/milc_qcd/blob/master/generic/ranstuff.c:27.
    Args:
        seed: int, seed for the random number generator
        index: int, selects algorithm / which random number generator / which multiplier
    """
    def __init__(self, seed, index):
        self.r = np.zeros(7, dtype=object)
        seed = (69607+8*index)*seed+12345
        self.r[0] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[1] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[2] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[3] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[4] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[5] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.r[6] = (seed>>8) & 0xffffff
        seed = (69607+8*index)*seed+12345
        self.ic_state = seed
        self.multiplier = 100005 + 8*index
        self.addend = 12345
        self.scale = 1.0/(0x1000000)

    def myrand(self):
        """
        Generates random numbers uniformly between zero and one.
        Returns:
            float, the random number
        """
        t = ( ((self.r[5] >> 7) | (self.r[6] << 17)) ^ ((self.r[4] >> 1) | (self.r[5] << 23)) ) & 0xffffff
        self.r[6] = self.r[5]
        self.r[5] = self.r[4]
        self.r[4] = self.r[3]
        self.r[3] = self.r[2]
        self.r[2] = self.r[1]
        self.r[1] = self.r[0]
        self.r[0] = t
        s = self.ic_state * self.multiplier + self.addend
        self.ic_state = s
        return self.scale*(t ^ ((s>>8)&0xffffff))

    def uniform(self, low=0.0, high=1.0, size=None):
        """
        Draw samples from a uniform distribution between "low" and "high".
        Args:
            low: float, the lower boundary of the output interval. Default is 0.0.
            high: float, the upper boundary of the output interval. Default is 1.0.
            size: int, the number of draws to make
        Returns:
            out: ndarray or scalar, the samples drawn
        """
        if high <= low:
            raise ValueError("Please specify high > low.")
        if size is None:
            size = 1
        out = np.array([self.myrand() for _ in np.arange(size)])
        out = out * (high - low) + low
        if size == 1:
            return out.item()
        return out


if __name__ == '__main__':
    main()