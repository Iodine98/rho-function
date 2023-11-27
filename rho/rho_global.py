from tqdm.auto import tqdm
RhoTuples = list[tuple[int, int]]

class RhoGlobal:

    RHO_MAP: dict[int, list[tuple[int, int]]] = {}
    PRIME_SET = set([2])

    def __init__(self, dynamic: bool = True, enable_tqdm: bool = False, trace: bool = False) -> None:
        self.dynamic = dynamic
        self.enable_tqdm = enable_tqdm
        self.trace = trace
    
    def break_down(self, n: int, base: int) -> tuple[int, list[tuple[int, int]]]:
        exponent = 0
        while n % base == 0:
            n = n // base
            exponent += 1
        if exponent > 0:
            return n, [(base, exponent)]
        return n, []
    
    def is_prime(self, k: int, add_to_prime_set: bool = False):
        for prime in self.PRIME_SET:
            if k % prime == 0 and prime != k:
                return False
        if add_to_prime_set:
            self.PRIME_SET.add(k)
        return True
    
    def rho_union(self, *args: RhoTuples, enable_tqdm: bool = False) -> RhoTuples:
        union_list = [*args[0]]
        progress_bar = tqdm(desc="list", total=len(args[1:]), disable=not enable_tqdm)
        for rho_right in args[1:]:
            for pair in rho_right:
                base, exp = pair
                indices = [i for i,x in enumerate(union_list) if base == x[0]]
                index = indices[0] if indices else None
                if index is not None:
                    union_list[index] = (union_list[index][0], union_list[index][1] + exp)
                    if union_list[index][1] == 0:
                        union_list.pop(index)
                else:
                    union_list.append(pair)
            progress_bar.update()
        progress_bar.close()
        return union_list
    
    def rho_power(self, rho_tuples: RhoTuples, power: int) -> RhoTuples:
        return [(base, exp * power) for base, exp in rho_tuples]
    

def count_decimals(number: float) -> int:
    """
    Count the number of decimals in a number.
    """
    number_as_str: str = str(number)
    if '.' in number_as_str:
        return len(number_as_str.split('.')[1])
    return 0