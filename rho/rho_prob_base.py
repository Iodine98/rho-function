from tqdm.auto import tqdm
import numpy as np
from rho.rho_global import RhoGlobal, RhoTuples, count_decimals


class RhoProbBase(RhoGlobal):
    def __init__(self, dynamic: bool = True, enable_tqdm: bool = False, trace: bool = False, threshold=np.int64(10**6)) -> None:
        super().__init__(dynamic, enable_tqdm, trace)
        self.threshold = threshold

    def __call__(self, n_input: np.int64, base=np.int64(2), inc=np.int64(1)) -> tuple[RhoTuples, np.float64]:
        return self.rho_prob_base(n_input, base, inc, self.enable_tqdm)

    def calc_prob(self, pair: tuple[np.int64, np.int64] | None) -> np.float64:
        """
        Calculate the probability of a pair.
        """
        return min(np.float64(1.00), (self.threshold / pair[0])) if pair else np.float64(1.00)

    # Define the Rho function
    def rho_prob_base(
        self, n_input: np.int64, base: np.int64, inc: np.int64, enable_tqdm: bool
    ) -> tuple[list[tuple[np.int64, np.int64]], np.float64]:
        """ "
        Rho base algorithm with probabilities for the likelihood of primeness.
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
            return pairs, self.calc_prob(pairs[-1] if pairs else None)

        while base < n_input:
            if base in self.PRIME_SET or self.is_prime(base, add_to_prime_set=True):
                if number == 1 or base > root or base > self.threshold:
                    break
                if number < n_input:
                    rho_bases, _ = self.rho_prob_base(number, base, inc, self.enable_tqdm and self.trace)
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
        return pairs, self.calc_prob(pairs[-1] if pairs else None)

    def rho_rational_number(self, f_input: float, base=np.int64(2), inc=np.int64(1)) -> tuple[RhoTuples, np.float64]:
        no_of_decimals = count_decimals(f_input)
        numerator = np.int64(round(f_input * 10**no_of_decimals))
        denominator = 10**no_of_decimals
        rho_numerator, p_numerator = self.rho_prob_base(numerator, base, inc, self.enable_tqdm)
        rho_denominator, _ = self.rho_prob_base(denominator, base, inc, self.enable_tqdm)
        rho = self.rho_union(rho_numerator, self.rho_power(rho_denominator, -1), enable_tqdm=False)
        return rho, p_numerator
