from tqdm.auto import tqdm
from rho.rho_global import RhoGlobal, RhoTuples, count_decimals
from decimal import Decimal


class RhoProbBase(RhoGlobal):
    def __init__(self, dynamic: bool = True, enable_tqdm: bool = False, trace: bool = False, threshold=10**6) -> None:
        super().__init__(dynamic, enable_tqdm, trace)
        self.threshold = threshold

    def __call__(self, n_input: int, base: int = 2, inc=1) -> tuple[RhoTuples, float]:
        return self.rho_prob_base(n_input, base, inc, self.enable_tqdm)

    def calc_prob(self, pair: tuple[int, int] | None) -> float:
        """
        Calculate the probability of a pair.
        """
        return min(1.00, (self.threshold / pair[0])) if pair else 1.00

    # Define the Rho function
    def rho_prob_base(self, n_input: int, base: int, inc: int, enable_tqdm: bool) -> tuple[list[tuple[int, int]], float]:
        """"
        Rho base algorithm with probabilities for the likelihood of primeness.
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
            return pairs, self.calc_prob(pairs[-1] if pairs else None)

        while base < n_input:
            if base in self.PRIME_SET or self.is_prime(base, add_to_prime_set=True):
                if number == 1 or base > root or base > self.threshold:
                    break
                if number < n_input:
                    rho_bases, _ = self.rho_prob_base(number, base, inc, self.enable_tqdm and self.trace)
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
        return pairs, self.calc_prob(pairs[-1] if pairs else None)

    def rho_rational_number(self, f_input: float, base: int = 2, inc=1) -> tuple[RhoTuples, float]:
        no_of_decimals = count_decimals(f_input)
        numerator = int(round(f_input * 10**no_of_decimals))
        denominator = 10**no_of_decimals
        rho_numerator, p_numerator = self.rho_prob_base(numerator, base, inc, self.enable_tqdm)
        rho_denominator, _ = self.rho_prob_base(denominator, base, inc, self.enable_tqdm)
        rho = self.rho_union(rho_numerator, self.rho_power(rho_denominator, -1), enable_tqdm=False)
        return rho, p_numerator