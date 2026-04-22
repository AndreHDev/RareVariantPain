
import os
import functools
import json
import hashlib

CACHE_FILE = "APIs/cache.txt"

def txt_cache(cache_file = CACHE_FILE):

    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, use_cache = True, **kwargs):

            if not use_cache:
                return func(*args, **kwargs)
            
            # Create deterministic cache key to recall a function that has already been called
            key_data = {
                "func": func.__name__,
                "args": args,
                "kwargs": kwargs
            }
            key_string = json.dumps(key_data, sort_keys=True)
            key_hash = hashlib.sha256(key_string.encode()).hexdigest()

            # Check cache
            if os.path.exists(cache_file):
                with open(cache_file, "r", encoding="utf-8") as f:
                    for line in f:
                        cached_key, cached_value = line.rstrip("\n").split("\t", 1)
                        if cached_key == key_hash:
                            print(f"{func.__name__}: Loaded from cache")
                            return json.loads(cached_value)

            # If not cached call function
            result = func(*args, **kwargs)

            # And append to cache
            with open(cache_file, "a", encoding="utf-8") as f:
                f.write(f"{key_hash}\t{json.dumps(result)}\n")

            print(f"{func.__name__}: Fetched from API and cached")
            return result

        return wrapper
    return decorator