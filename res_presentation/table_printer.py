def _get_rows_from_row_with_lists(row):
    def get_max_size(lst):
        max_size = 0
        for item in lst:
            if isinstance(item, list):
                size = len(item)
                if size > max_size:
                    max_size = size
        return max_size
    
    def transpose_matrix(matrix):
        transposed = [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]
        return transposed

    max_size = get_max_size(row)
    rows_with_expanded_columns = []
    for item in row:
        if isinstance(item, list):
            diff = max_size - len(item)
            padded_list = item + [None] * diff
            rows_with_expanded_columns.append(padded_list)
        else:
            rows_with_expanded_columns.append([item] + [None] * (max_size - 1))

    rows = transpose_matrix(rows_with_expanded_columns)
    return rows


def _row_items_to_strings(multi_row):
    def item_to_string(item):
        if item is None:
            return ''
        if isinstance(item, (int, float)):
            return "{:.3f}".format(item)
        return item

    return [
        [ item_to_string(item) for item in row ]
        for row in multi_row
    ]

def _normalize(table):
    multi_rows = [
        _get_rows_from_row_with_lists(multi_row)
        for multi_row in table
    ]
    
    normalized_stringified = [
        _row_items_to_strings(multi_row)
        for multi_row in multi_rows
    ]
    
    return normalized_stringified

def _get_maxes_by_columns(table):
    def get_max_for_sub_table(sub_table, column):
        return max([len(row[column]) for row in sub_table])
    
    def get_max_for_column(table, column):
        return max([
            get_max_for_sub_table(sub_table, column)
            for sub_table in table
        ])

    columns = len(table[0][0])

    return [
        get_max_for_column(table, i)
        for i in range(columns)
    ]

def _print_sub_table(sub_table, columns_widths):
    def get_part(part, desired_width):
        return part + ' ' * (desired_width - len(part))
    
    def print_line(row):
        parts = [
            get_part(row[col], columns_widths[col])
            for col in range(len(columns_widths))
        ]

        print(f"| {' | '.join(parts)} |")

    for row in sub_table:
        print_line(row)

def _print_horizontal_line(columns_widths):
    print(f"+{'+'.join(['-' * (w + 2) for w in columns_widths])}+")

def print_table(table):
    normalized_table = _normalize(table)
    columns_widths = _get_maxes_by_columns(normalized_table)
    
    for sub_table in normalized_table:
        _print_horizontal_line(columns_widths)
        _print_sub_table(sub_table, columns_widths)
            
    _print_horizontal_line(columns_widths)
    