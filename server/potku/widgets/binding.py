# coding=utf-8
"""
Created on 18.03.2020

Potku is a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2020 TODO

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (file named 'LICENCE').
"""
__author__ = "Juhani Sundell"
__version__ = "2.0"

import os
import abc
import json
import time

from widgets.scientific_spinbox import ScientificSpinBox
from widgets.isotope_selection import IsotopeSelectionWidget

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import QTime


def from_qtime(qtime: QTime) -> int:
    """Converts QTime object to seconds.
    """
    return qtime.hour() * 60 * 60 + qtime.minute() * 60 + qtime.second()


def to_qtime(seconds: int) -> QTime:
    """Converts seconds to QTime.
    """
    t = QTime(0, 0, 0).addSecs(seconds)
    return t


def unix_time_to_label(instance, attr: str, unix_time: float):
    """Sets the text of a label to a formatted string representations
    of the given unix time. Unix time is also stored as an attribute of
    the label so it can be later retrieved.
    """
    label = getattr(instance, attr)
    label.unix_time = unix_time
    label.setText(time.strftime("%c %z %Z", time.localtime(unix_time)))


def unix_time_from_label(instance, attr: str) -> float:
    label = getattr(instance, attr)
    try:
        return label.unix_time
    except AttributeError:
        return 0.0


def get_items(list_widget: QtWidgets.QListWidget):
    """Returns the items in a QListWidget.
    """
    return [
        list_widget.item(i).data(QtCore.Qt.UserRole)
        for i in range(list_widget.count())
    ]


def set_items(list_widget: QtWidgets.QListWidget, ls):
    """Sets the items in a QListWidget.
    """
    list_widget.clear()
    for li in ls:
        item = QtWidgets.QListWidgetItem(str(li))
        item.setData(QtCore.Qt.UserRole, li)
        list_widget.addItem(item)


def set_selected_combobox_item(combobox: QtWidgets.QComboBox, data):
    """Sets the current item in a combobox depending on the given
    data.
    """
    idx = combobox.findData(data, QtCore.Qt.UserRole, QtCore.Qt.MatchExactly)
    if idx != -1:
        combobox.setCurrentIndex(idx)


def get_btn_group_value(button_group: QtWidgets.QButtonGroup):
    """Retuns the value of data_item attribute of the selected button.
    """
    try:
        return button_group.checkedButton().data_item
    except AttributeError:
        return None


def set_btn_group_value(button_group: QtWidgets.QButtonGroup, value):
    """Checks the button whose data_item matches the value.
    """
    for btn in button_group.buttons():
        if btn.data_item == value:
            btn.setChecked(True)
            break


# Collections of default getter and setter methods for various QObjects.
# Keys are the types of the QObjects and values are methods.
_DEFAULT_GETTERS = {
    QtWidgets.QTimeEdit: lambda qobj: from_qtime(qobj.time()),
    QtWidgets.QLineEdit: lambda qobj: qobj.text(),
    QtWidgets.QComboBox: lambda qobj: qobj.currentData(QtCore.Qt.UserRole),
    QtWidgets.QTextEdit: lambda qobj: qobj.toPlainText(),
    QtWidgets.QCheckBox: lambda qobj: qobj.isChecked(),
    QtWidgets.QLabel: lambda qobj: qobj.text(),
    QtWidgets.QPlainTextEdit: lambda qobj: qobj.toPlainText(),
    QtWidgets.QPushButton: lambda qobj: qobj.text(),
    ScientificSpinBox: lambda qobj: qobj.get_value(),
    IsotopeSelectionWidget: lambda qobj: qobj.get_element(),
    QtWidgets.QListWidget: get_items,
    QtWidgets.QGroupBox: lambda qobj: qobj.title(),
    QtWidgets.QButtonGroup: get_btn_group_value
}

_DEFAULT_SETTERS = {
    QtWidgets.QTimeEdit: lambda qobj, sec: qobj.setTime(to_qtime(sec)),
    QtWidgets.QLineEdit: lambda qobj, txt: qobj.setText(txt),
    QtWidgets.QComboBox: set_selected_combobox_item,
    QtWidgets.QTextEdit: lambda qobj, txt: qobj.setText(txt),
    QtWidgets.QCheckBox: lambda qobj, b: qobj.setChecked(b),
    QtWidgets.QLabel: lambda qobj, txt: qobj.setText(str(txt)),
    QtWidgets.QPlainTextEdit: lambda qobj, txt: qobj.setPlainText(txt),
    QtWidgets.QPushButton: lambda qobj, txt: qobj.setText(txt),
    ScientificSpinBox: lambda qobj, value: qobj.set_value(value),
    IsotopeSelectionWidget: lambda qobj, elem: qobj.set_element(elem),
    QtWidgets.QListWidget: set_items,
    QtWidgets.QGroupBox: lambda qobj, txt: qobj.setTitle(txt),
    QtWidgets.QButtonGroup: set_btn_group_value
}


def _fget(instance, qobj_name):
    """Returns the value of a QObject.

    Args:
        instance: object that holds a reference to a QObject.
        qobj_name: name of the reference to a QObject.

    Return:
        value of the QObject.
    """
    qobj = getattr(instance, qobj_name)
    getter = _DEFAULT_GETTERS.get(type(qobj), lambda obj: obj.value())
    return getter(qobj)


def _fset(instance, qobj_name, value):
    """Sets the value of a QObject.

    Args:
        instance: object that holds a reference to a QObject.
        qobj_name: name of the reference to a QObject.
        value: new value for the QObject.
    """
    qobj = getattr(instance, qobj_name)
    setter = _DEFAULT_SETTERS.get(
        type(qobj), lambda obj, val: obj.setValue(val))
    setter(qobj, value)


class PropertyBindingWidget(abc.ABC):
    """Base class for a widget that contains bindable properties.
    """
    def _get_properties(self, only_tracking_properties=False):
        """Returns a generator of the names of all the properties the widget
        has.

        Args:
            only_tracking_properties: whether only TrackingProperties are
                yielded.
        """
        if only_tracking_properties:
            def cond(prop):
                return isinstance(getattr(type(self), prop), TrackingProperty)
        else:
            def cond(prop):
                return isinstance(getattr(type(self), prop), property)
        return (
            d for d in dir(self)
            if hasattr(type(self), d) and
            cond(d)
        )

    def set_properties(self, **kwargs):
        """Sets property values.

        Args:
            kwargs: property names and values as keyword arguments.
        """
        for p in self._get_properties():
            if p in kwargs:
                try:
                    setattr(self, p, kwargs[p])
                except (AttributeError, TypeError):
                    pass

    def get_properties(self):
        """Returns property names and their values as a dictionary.

        Return:
            a dictionary.
        """
        return {
            p: getattr(self, p)
            for p in self._get_properties()
        }


class PropertyTrackingWidget(PropertyBindingWidget, abc.ABC):
    """Widget that stores the original values of its properties, and is
    able to check if the values have changed.
    """
    # TODO possibly add methods for resetting and clearing original property
    #  values.

    @abc.abstractmethod
    def get_original_property_values(self):
        """Returns a dictionary of property names and their values.

        The purpose of this function is to provide a dictionary for a
        TrackingProperty to store values.
        """
        # TODO make this an underscore method
        pass

    def are_values_changed(self):
        """Checks if current property values differ from original ones.

        Return:
            boolean.
        """
        return any(getattr(type(self), p).is_value_changed(self)
                   for p in self._get_properties(only_tracking_properties=True))


class PropertySavingWidget(PropertyBindingWidget, abc.ABC):
    """Property widget that saves the current state of its properties
    to file. Saving is done automatically when the widget closes, unless
    the inheriting widget overrides the closeEvent function.

    Loading is done by calling the load_properties_from_file function.

    Property values must be JSON encodable.
    """
    # TODO maybe do loading in showEvent function

    @abc.abstractmethod
    def get_property_file_path(self):
        """Returns Path object to the file that is used to save and load
        properties.
        """
        pass

    def closeEvent(self, event):
        """Overrides QWidgets closeEvent. Saves the properties to default
        file.
        """
        self.save_properties_to_file(file_path=self.get_property_file_path())
        event.accept()

    def save_properties_to_file(self, file_path=None):
        """Saves properties to a file.

        Args:
            file_path: Path object to a file that is used for saving. If None,
                widget's default property file location is used.
        """
        if file_path is None:
            file_path = self.get_property_file_path()

        os.makedirs(file_path.parent, exist_ok=True)
        params = self.get_properties()
        try:
            with open(file_path, "w") as file:
                json.dump(params, file, indent=4)
        except (PermissionError, IsADirectoryError):
            pass

    def load_properties_from_file(self, file_path=None):
        """Loads properties from a file.

        Args:
            file_path: Path object to a file that is used for loading. If None,
                widget's default property file location is used.
        """
        if file_path is None:
            file_path = self.get_property_file_path()
        try:
            with open(file_path) as file:
                params = json.load(file)
        except (OSError, json.JSONDecodeError, UnicodeDecodeError,
                IsADirectoryError):
            return

        self.set_properties(**params)


class TrackingProperty(property):
    """TrackingProperty can be used in a PropertyTrackingWidget to detect
    if the property value has changed from the first time the property was set.

    If TrackingProperty is used in another kind of object, it behaves like
    a normal property.
    """
    # Note: the reason why original value is stored in the widget instead of
    # the property, is because the property is a class attribute and therefore
    # common for all instances of that class.
    # Alternatively, the property could store the original values and a weak
    # reference to every instance in a dictionary.

    def __init__(self, *args, attr_name=None, **kwargs):
        """Initializes a TrackingProperty.

        Args:
            args: arguments passed down to property constructor.
            kwargs: keyword arguments passed down to property constructor.
            attr_name: name of the attribute that the property is bound to. If
                value is None, the value of the property is not tracked and
                TrackingProperty behaves like a normal property.
        """
        super().__init__(*args, **kwargs)
        self.__prop_name = attr_name

    def __set__(self, instance, value):
        """Sets the value of the property.

        Args:
            instance: object that holds a reference to the property.
            value: new value of the property.
        """
        # TODO consider using a pyQtSignal when setting values of QObjects
        #   to ensure that GUI is updated in the main thread.
        self.fset(instance, value)
        if self.__prop_name is not None and \
                isinstance(instance, PropertyTrackingWidget):
            # Only PropertyTrackingWidget can store original values.

            orig_props = instance.get_original_property_values()
            if self.__prop_name not in orig_props:
                # If the property has not yet been stored, store it now.
                orig_props[self.__prop_name] = self.fget(instance)

    def is_value_changed(self, instance):
        """Checks if the current value of the property differs from the original
        one.

        Args:
            instance: object that has a reference to the property. If the
                object is not an instance of PropertyTrackingWidget,
                this returns False.

        Return:
            boolean.
        """
        if self.__prop_name is not None and \
                isinstance(instance, PropertyTrackingWidget):
            cur_value = self.fget(instance)
            orig_props = instance.get_original_property_values()
            orig_value = orig_props.get(self.__prop_name, cur_value)
            return cur_value != orig_value
        return False


def bind(attr_name, fget=None, fset=None, twoway=True,
         track_change=False):
    """Returns a property that is bound to an attribute.

    Mostly used to bind a QObject value to a property, but other attributes
    can be used as well if custom getter and setter are defined.

    Args:
        attr_name: name of an attribute that the property will be bound to
        fget: getter function for the property. If None, a default getter
            function is used depending on the type of QObject that the attr_name
            references. Functions should have the signature:
                f(instance, attr_name)
        fset: setter function for the property. If None, a default setter
            function is used depending on the type of QObject that the attr_name
            references. Functions should have the signature:
                f(instance, attr_name, value)
        twoway: whether binding is two-way (setting the property also sets
            the QObject value) or one-way (only getter is defined for the
            property).
        track_change: whether the property tracks changes in its value or not.

    Return:
        TrackingProperty.
    """
    def getter(instance):
        if fget is None:
            return _fget(instance, attr_name)
        return fget(instance, attr_name)

    if not twoway:
        return TrackingProperty(getter)

    def setter(instance, value):
        try:
            if fset is None:
                _fset(instance, attr_name, value)
            else:
                fset(instance, attr_name, value)
        except TypeError:
            # value is wrong type, nothing to do
            pass

    if not track_change:
        # If changes in property names do not need to be tracked, pass the
        # name as None for the TrackingProperty
        prop_name_ = None
    else:
        prop_name_ = attr_name

    return TrackingProperty(getter, setter, attr_name=prop_name_)


def multi_bind(attrs, fget=None, fset=None, twoway=True, track_change=False):
    # TODO refactor this with bind
    # TODO enable twoway binding with combobox
    def getter(instance):
        if fget is None:
            return tuple(
                _fget(instance, attr)
                for attr in attrs
            )
        else:
            return fget(instance, attrs)

    if not twoway:
        return TrackingProperty(getter)

    def setter(instance, values):
        if fset is None:
            for attr, value in zip(attrs, values):
                try:
                    _fset(instance, attr, value)
                except TypeError:
                    pass
        else:
            fset(instance, attrs, values)

    if track_change:
        attrs_ = tuple(attrs)
    else:
        attrs_ = None

    return TrackingProperty(getter, setter, attr_name=attrs_)


