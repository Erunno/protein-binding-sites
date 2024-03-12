import json
import os
import pprint
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import config.config as config
from typing import TypeVar

class FileCache:
    def __init__(self, file_key, data_folder):
        self.file = os.path.join(data_folder, f'{file_key}.json')
        self.cache = None

        self.load()
    
    def load(self):
        try: 
            with open(self.file, 'r') as file:
                self.cache = json.load(file)
        except:
            self.cache = {}

    def get_value_for(self, method_name, *args, **kwargs):
        if method_name not in self.cache:
            return False, None
        
        cache_record = self.cache[method_name]
        key_to_cache = self.__get_key_to_cache(*args, **kwargs)

        if key_to_cache not in cache_record:
            return False, None
        
        return True, cache_record[key_to_cache]

    def save(self, value, method_name, *args, **kwargs):
        if method_name not in self.cache:
            self.cache[method_name] = {}
        
        cache_record = self.cache[method_name]
        key_to_cache = self.__get_key_to_cache(*args, **kwargs)

        cache_record[key_to_cache] = value

    def flush(self):
        result = pprint.pformat(self.cache, compact=True, width=150).replace("'",'"')
        with open(self.file, 'w') as file:
            file.write(result)
        
    def free_memory(self):
        self.cache = None

    def __get_key_to_cache(self, *args, **kwargs):
        args_str = ', '.join(repr(arg) for arg in args)
        kwargs_str = ', '.join(
            sorted([f"{key}={repr(value)}" for key, value in kwargs.items()])
        )
        return f"{args_str}, {kwargs_str}"
                
T = TypeVar('T')
def use_cache(object: T, data_folder=config.cache_folder) -> T:
    cache = FileCache(file_key=object.__cache_key__, data_folder=data_folder)

    class CachedCalls:
        def __enter__(self):
            return self
        
        def __exit__(self, exc_type, exc_value, traceback):
            cache.flush()
            cache.free_memory()

        def __getattribute__(self, name):
            attribute = object.__getattribute__(name)
            cashable_attributes = object.__getattribute__('__cashable__')
            
            if name not in cashable_attributes:
                return attribute 

            method = attribute

            def cached_call_wrapper(*args, **kwargs):
                is_in_cache, value = cache.get_value_for(
                    name, *args, **kwargs
                )

                if is_in_cache:
                    return value

                result = method(*args, **kwargs)
                cache.save(result, name, *args, **kwargs)

                return result

            return cached_call_wrapper

    return CachedCalls()
