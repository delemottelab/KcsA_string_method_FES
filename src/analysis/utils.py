
def natural_sort(my_list):
    """
    Takes as input a list l of strings and sorts it with natural order.
      Parameters
      ----------
      l: list of strings.
      Returns
      -------
      l sorted
    """
    from re import split

    assert isinstance(my_list, list), "l is not a list!"
    for i in my_list:
        assert isinstance(i, str), "List contains non-string elements."

    def convert(text):

        if text.isdigit():
            convert = int(text)
        else:
            convert = text.lower()
        return convert

    def alphanum_key(key):

        return [convert(c) for c in split("([0-9]+)", key)]

    return sorted(my_list, key=alphanum_key)
