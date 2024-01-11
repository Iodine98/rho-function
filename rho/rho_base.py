from typing import Any
from decimal import Decimal

from tqdm.auto import tqdm

from rho.rho_global import RhoGlobal, RhoTuples, count_decimals


class RhoBase(RhoGlobal):
    """
    Rho base algorithm.
    """

    def __call__(self, n_input: int, base: int = 2, inc=1) -> Any:
        return self.rho_base(n_input, base, inc, self.enable_tqdm)

    def rho_base(self, n_input: int, base: int, inc: int, enable_tqdm: bool) -> RhoTuples:
        """
        Rho base algorithm.
        """
        pairs: list[tuple[int, int]] = self.RHO_MAP.get(n_input, []) if self.dynamic else []
        number = n_input
        range_bar = tqdm(
            total=n_input,
            initial=base,
            disable=not enable_tqdm,
            unit_scale=True,
        )
        root = Decimal(n_input).sqrt()

        if pairs:
            return pairs

        while base < n_input:
            if base in self.PRIME_SET or self.is_prime(base, add_to_prime_set=True):
                if number == 1 or base > root:
                    break
                if number < n_input:
                    rho_bases = self.rho_base(number, base, inc, self.enable_tqdm and self.trace)
                    pairs.extend(rho_bases)
                    number = 1
                    continue
                number, pair = self.break_down(number, base)
                pairs.extend(pair)
            base += inc
            range_bar.update(inc)
        range_bar.update(max(0, n_input - base))
        if number > 1:
            pairs.append((number, 1))
        range_bar.close()
        return pairs

    def rho_rational_number(self, f_input: float, base: int = 2, inc=1) -> RhoTuples:
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
        numerator = int(round(f_input * 10**no_of_decimals))
        denominator = 10**no_of_decimals
        rho_numerator = self.rho_base(numerator, base, inc, self.enable_tqdm)
        rho_denominator = self.rho_base(denominator, base, inc, self.enable_tqdm)
        rho = self.rho_union(rho_numerator, self.rho_power(rho_denominator, -1), enable_tqdm=False)
        return rho
