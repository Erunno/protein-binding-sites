from typing import TypeVar

T = TypeVar('T')

def use_cache(object: T) -> T:
    cache = {}

    class CachedCalls:
        def __enter__(self):
            return self
        
        def __exit__(self, exc_type, exc_value, traceback):
            print ('cache was', cache)
            return  

        def __getattribute__(self, name):
            attribute = object.__getattribute__(name)
            cashable = object.__getattribute__('__cashable__')
            
            if name not in cashable:
                return attribute 

            method = attribute

            def cached_call_wrapper(*args, **kwargs):
                def to_key(*args, **kwargs):
                    args_str = ', '.join(repr(arg) for arg in args)
                    kwargs_str = ', '.join(f"{key}={repr(value)}" for key, value in kwargs.items())

                    result_str = f"{name}({args_str}, {kwargs_str})"
                    return result_str
                
                key_to_cache = to_key(*args, **kwargs)
                print ('key', key_to_cache)

                if key_to_cache in cache:
                    return cache[key_to_cache]

                result = method(*args, **kwargs)


                cache[key_to_cache] = result
                return result
            

            return cached_call_wrapper

    return CachedCalls()

# Example usage:

class BaseClass:
    def __init__(self):
        self.__cache_key__ = 'key_to_cache'
        self.__cashable__ = [
            self.expensive_operation.__name__
        ]

    def expensive_operation(self, *args, **kwargs):
        print(f"Performing expensive operation... for {args}, {kwargs}")
        return f"Result of expensive operation ... for {args}, {kwargs}"
    

base_class = BaseClass()

with use_cache(base_class) as cached_base:

    res1 = cached_base.expensive_operation(1)
    print()
    res2 = cached_base.expensive_operation(1, test='test')
    print()
    res3 = cached_base.expensive_operation(1, test='test')
    print()

    print ('res 1:', res1)
    print ('res 2:', res2)
    print ('res 3:', res3)
