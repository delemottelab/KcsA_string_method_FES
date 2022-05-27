import numpy as np


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


def find_nearest_point(data, point):
    data_ = data - point
    data_ = np.sum(data_ * data_, axis=1)
    nearest = data_.argmin()
    distance = np.sqrt(data_.min())
    return nearest, distance


def jupyter_lab_notification(
    time_in_seconds=180, sound_path="/usr/share/sounds/gnome/default/alerts/sonar.ogg"
):
    from jupyter_helpers.notifications import Notifications

    Notifications(
        success_audio=sound_path, time_threshold=time_in_seconds, integration="GNOME"
    )
    return


def jupyter_lab_error(sound_path="/usr/share/sounds/gnome/default/alerts/bark.ogg"):
    from jupyter_helpers.notifications import Notifications

    Notifications(failure_audio=sound_path, integration="GNOME")

    return
