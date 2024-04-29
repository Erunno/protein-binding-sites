import json
import os

def load_from(folder, limit=None):
    results = []
    file_names = [fname for fname in os.listdir(folder) if fname.endswith('.json')]

    total = len(file_names) if not limit else limit

    for i, filename in enumerate(file_names[:limit]):
        if i % 100:
            print(f'\r[{"=" * ((i + 1) * 30 // total)}{" " * ((total - i - 1) * 30 // total)}]({i + 1:5}/{total:5})', end='')
        
        file_path = os.path.join(folder, filename)

        with open(file_path, 'r') as file:
            try:
                result = json.load(file)
            except:
                print (f'ERROR: cannot read file: {file_path}')
                exit(1)

            result['file_name'] = filename
            results.append(result)

    print(' ' * 50)
    return results