from typing import Any
import numpy as np
from tqdm.auto import tqdm

from rho.rho_global import RhoGlobal, RhoTuples, count_decimals


class RhoBase(RhoGlobal):
    """
    Rho base algorithm.
    """

    def __call__(self, n_input: np.int64, base=np.int64(2), inc=np.int64(1)) -> Any:
        return self.rho_base(n_input, base, inc, self.enable_tqdm)

    def rho_base(self, n_input: np.int64, base: np.int64, inc: np.int64, enable_tqdm: bool) -> RhoTuples:
        """
        Rho base algorithm.
        """
        pairs: list[tuple[np.int64, np.int64]] = self.RHO_MAP.get(n_input, []) if self.dynamic else []
        number = n_input
        range_bar = tqdm(
            total=int(n_input),
            initial=int(base),
            disable=not enable_tqdm,
            unit_scale=True,
        )
        root = np.sqrt(n_input)

        if pairs:
            return pairs

        while base < n_input:
            if base in self.PRIME_SET or self.is_prime(base, add_to_prime_set=True):
                if number == 1 or base > root:
                    break
                if number < n_input:
                    rho_bases = self.rho_base(number, base, inc, self.enable_tqdm and self.trace)
                    pairs.extend(rho_bases)
                    number = np.int64(1)
                    continue
                number, pair = self.break_down(number, base)
                pairs.extend(pair)
            base += inc
            range_bar.update(int(inc))
        range_bar.update(int(max(0, n_input - base)))
        if number > 1:
            pairs.append((number, np.int64(1)))
        range_bar.close()
        if self.dynamic and n_input not in self.RHO_MAP:
            self.RHO_MAP[n_input] = pairs
        return pairs

    def rho_rational_number(self, f_input: float, base: np.int64 = np.int64(2), inc=np.int64(1)) -> RhoTuples:
        """
        Calculates the Rho value and the probability for a given rational number.

        Args:
            f_input (float): The input rational number.
            base (int, optional): The base for the Rho calculation. Defaults to 2.
            inc (int, optional): The increment value for the Rho calculation. Defaults to 1.

        Returns:
            tuple[RhoTuples, float]: A tuple containing the Rho value and the probability.
        """
        no_of_decimals = count_decimals(f_input)
        numerator = np.int64(round(f_input * 10**no_of_decimals))
        denominator = 10**no_of_decimals
        rho_numerator = self.rho_base(numerator, base, inc, self.enable_tqdm)
        rho_denominator = self.rho_base(denominator, base, inc, self.enable_tqdm)
        rho = self.rho_union(rho_numerator, self.rho_power(rho_denominator, -1), enable_tqdm=False)
        return rho
