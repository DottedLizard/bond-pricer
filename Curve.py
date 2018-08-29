import datetime
import re
from Date import Date
import collections
import math
from Data_Utils import DataUtils


class Curve(object):
    """Base class for curve objects"""

    _valid_interpolation_methods = ["ZL"]
    _valid_frequencies = ["CONTINUOUS", "CONT", "ANNUAL", "A", "1", "SEMI-ANNUAL", "S", "2", "QUARTERLY", "Q", "4",
                          "MONTHLY", "M", "12", "DAILY", "D", "365", "SIMPLE", "SMP"]
    _valid_day_count_conv = ['1/1', 'ACT/ACT', 'ACT/365', 'ACT/360', '30/360', '30E/360']

    def __init__(self, name, base_date, dates, values, interpolation="ZL", frequency="CONT", day_count_conv="ACT/365"):

        self._name = name
        self._base_date = base_date
        self._dates = dates
        self._normalized_dates = []
        self._values = values
        self._interpolation = interpolation
        self._frequency = frequency.upper()
        self._day_count_conv = day_count_conv.upper()
        self._curve = collections.OrderedDict()

        self._validate_inputs()
        self._construct_curve()

    @property
    def name(self):
        return self._name

    @property
    def base_date(self):
        return self._base_date

    @property
    def dates(self):
        return self._normalized_dates

    @property
    def values(self):
        return self._values

    @property
    def curve(self):
        return self._curve

    def _validate_inputs(self):

        DataUtils.validate_type(self._name, str, "Name")
        DataUtils.validate_type(self._base_date, datetime.date, "Base date")
        DataUtils.validate_type(self._dates, list, "Dates")
        DataUtils.validate_type(self._values, list, "Values")
        self._validate_dates()
        DataUtils.validate_same_length([self.dates, self.values], ["Date list", "Value list"])
        DataUtils.validate_in_list(self._interpolation, Curve._valid_interpolation_methods, "Interpolation")
        DataUtils.validate_in_list(self._frequency, Curve._valid_frequencies, "Compounding frequency")
        DataUtils.validate_in_list(self._day_count_conv, Curve._valid_day_count_conv, "Day count convention")

    def _validate_dates(self):

        if not (self._all_datetime() or self._all_tenors()):
            raise TypeError("Dates in date list must be of type datetime.date or tenor strings")

        self._normalize_dates()

        if min(self.dates) < self.base_date:
            raise ValueError("Dates must begin on or after base date")

        if not all(x < y for x, y in zip(self.dates, self.dates[1:])):
            raise ValueError("Dates must be strictly increasing")

    def _normalize_dates(self):

        if self._all_datetime():
            self._normalized_dates = self._dates

        else:

            for tenor in self._dates:

                normalized_date = Date.add_tenor(self.base_date, tenor)
                self._normalized_dates.append(normalized_date)

    def _all_datetime(self):
        return all([(type(date) is datetime.date) for date in self._dates])

    def _all_tenors(self):
        return all([(re.match("^\d+[DMY]$", str(date)) is not None) for date in self._dates])

    def _construct_curve(self):

        for i in range(0, len(self.dates)):
            self._curve[self.dates[i]] = self.values[i]

    def spot_rate(self, date):

        DataUtils.validate_date_or_tenor(date)

        if date < self.base_date:
            raise ValueError("Date must be greater than or equal to base date")

        ref_date = self._calc_ref_date(date)

        if ref_date in self.dates:
            return self._curve[ref_date]
        else:
            return self._interpolated_curve_value(ref_date)

    def _calc_ref_date(self, date):

        if re.match("^\d+[DMY]$", str(date)) is not None:
            return Date.add_tenor(self.base_date, date)
        else:
            return date

    def _interpolated_curve_value(self, ref_date):

        if self._interpolation == "ZL":

            if ref_date <= min(self.dates):
                return self.values[0]
            elif ref_date >= max(self.dates):
                return self.values[len(self.values) - 1]

            else:

                start_idx = 0
                end_idx = 0

                for i in range(0, len(self.dates)):

                    if ref_date <= self.dates[i]:

                        start_idx = i - 1
                        end_idx = i
                        break

                total_period_days = (self.dates[end_idx] - self.dates[start_idx]).days
                elapsed_period_days = (ref_date - self.dates[start_idx]).days
                value_increment = self.values[end_idx] - self.values[start_idx]

                return self.values[start_idx] + ((elapsed_period_days / total_period_days) * value_increment)

        else:

            raise ValueError("Invalid interpolation method")

    def parallel_shift(self, shift):

        curve = YieldCurve(self._name, self._base_date, self._dates, self._values, self._interpolation, self._frequency,
                           self._day_count_conv)

        for key, value in curve._curve.items():
            curve._curve[key] = value + shift

        curve._values = list(curve._curve.values())

        return curve


class SwapCurve(Curve):
    """Class to create and manipulate swap curve objects"""

    def __init__(self, name, base_date, dates, values, interpolation="ZL", frequency="S",
                 day_count_conv="ACT/365"):

        super().__init__(name, base_date, dates, values, interpolation, frequency, day_count_conv)


class YieldCurve(Curve):
    """Class to create and manipulate yield curve objects"""

    def __init__(self, name, base_date, dates, values, interpolation="ZL", frequency="CONT",
                 day_count_conv="ACT/365"):

        super().__init__(name, base_date, dates, values, interpolation, frequency, day_count_conv)

        self._compounding_periods = {
            "ANNUAL": 1,
            "A": 1,
            "1": 1,
            "SEMI-ANNUAL": 2,
            "S": 2,
            "2": 2,
            "QUARTERLY": 4,
            "Q": 4,
            "4": 4,
            "MONTHLY": 12,
            "M": 12,
            "12": 12,
            "DAILY": 365,
            "D": 365,
            "365": 365
        }

    def df(self, date):

        DataUtils.validate_type(date, datetime.date, "Date")

        rate = self.spot_rate(date)
        year_frac = Date.year_frac(self.base_date, date, self._day_count_conv)

        return self.__calc_df(rate, year_frac)

    def fwd_df(self, date1, date2):

        DataUtils.validate_date_or_tenor(date1)
        DataUtils.validate_date_or_tenor(date2)

        ref_date1 = self._calc_ref_date(date1)
        ref_date2 = self._calc_ref_date(date2)

        fwd_rate = self.fwd_rate(ref_date1, ref_date2)

        year_frac = Date.year_frac(date1, date2, self._day_count_conv)

        return self.__calc_df(fwd_rate, year_frac)

    def __calc_df(self, rate, year_frac):

        if (self._frequency == "CONT") or (self._frequency == "CONTINUOUS"):
            return math.exp(-rate * year_frac)
        elif (self._frequency == "SMP") or (self._frequency == "SIMPLE"):
            return 1 / (1 + rate * year_frac)
        else:

            n = self._compounding_periods[self._frequency]

            return (1 + (rate / n)) ** (-year_frac * n)

    def convert_to_convention(self, name, frequency, day_count_conv):

        frequency = frequency.upper()
        day_count_conv = day_count_conv.upper()

        DataUtils.validate_type(name, str, "Name")
        DataUtils.validate_in_list(frequency, YieldCurve._valid_frequencies, "Compounding frequency")
        DataUtils.validate_in_list(day_count_conv, YieldCurve._valid_day_count_conv, "Day count convention")

        new_values = []

        for key, value in self.curve.items():

            accum_factor = 1 / self.df(key)
            year_frac = Date.year_frac(self.base_date, key, day_count_conv)

            if year_frac == 0:
                new_values.append(value)

            else:

                if (frequency == "CONT") or (frequency == "CONTINUOUS"):
                    new_values.append(math.log(accum_factor**(1 / year_frac)))

                elif (frequency == "SMP") or (frequency == "SIMPLE"):
                    new_values.append((accum_factor - 1) / year_frac)

                else:
                    new_values.append(self._compounding_periods[frequency] *
                                      (accum_factor**(1 / (year_frac * self._compounding_periods[frequency])) - 1))

        new_yield_curve = YieldCurve(name, self.base_date, self.dates, new_values, self._interpolation,
                                     frequency, day_count_conv)
        return new_yield_curve

    def fwd_rate(self, date1, date2):

        DataUtils.validate_date_or_tenor(date1)
        DataUtils.validate_date_or_tenor(date2)

        ref_date1 = self._calc_ref_date(date1)
        ref_date2 = self._calc_ref_date(date2)

        if ref_date1 < self.base_date:
            raise ValueError("Start date must be greater than or equal to base date")

        if ref_date2 <= ref_date1:
            raise ValueError("End date must be greater than or equal to start date")

        rate1 = self.spot_rate(date1)
        rate2 = self.spot_rate(date2)

        year_frac1 = Date.year_frac(self.base_date, ref_date1, self._day_count_conv)
        year_frac2 = Date.year_frac(self.base_date, ref_date2, self._day_count_conv)

        if (self._frequency == "CONT") or (self._frequency == "CONTINUOUS"):
            return ((rate2 * year_frac2) - (rate1 * year_frac1)) / (year_frac2 - year_frac1)
        elif (self._frequency == "SMP") or (self._frequency == "SIMPLE"):
            return ((1 + rate2 * year_frac2) / (1 + rate1 * year_frac1) - 1) / (year_frac2 - year_frac1)
        else:

            n = self._compounding_periods[self._frequency]
            compounded_fv_rate1 = (1 + (rate1 / n))**(year_frac1 * n)
            compounded_fv_rate2 = (1 + (rate2 / n))**(year_frac2 * n)

            if (year_frac2 - year_frac1) == 0.:
                return 0.

            annualizing_exp = 1 / ((year_frac2 - year_frac1) * n)

            return n * ((compounded_fv_rate2 / compounded_fv_rate1)**annualizing_exp - 1)

'''
base_date = datetime.date(2016, 10, 25)
dates = []
values = []

for i in range(0, 100):

    tenor = str(i) + "M"
    dates.append(Date.add_tenor(base_date, tenor))
    values.append(0.001 + i * 0.0004)

swap_curve = SwapCurve("test_curve", base_date, dates, values)
ref_date = datetime.date(2020, 12, 31)
print(swap_curve.spot_rate(ref_date))
print(swap_curve.name)
print(swap_curve.base_date)
print(swap_curve.dates)
print(swap_curve.values)
print(swap_curve.curve)

yield_curve = YieldCurve("yield_curve", base_date, dates, values, "ZL", "M", "Act/Act")
print("YIELD CURVE:")
print(yield_curve.spot_rate(ref_date))
print(yield_curve.name)
print(yield_curve.base_date)
print(yield_curve.dates)
print(yield_curve.values)
print(yield_curve.curve)

date1 = datetime.date(2016, 12, 22)
date2 = datetime.date(2017, 1, 22)
date3 = datetime.date(2017, 2, 28)
date4 = datetime.date(2017, 6, 12)
date5 = datetime.date(2020, 4, 3)
print(yield_curve.df(date1))
print(yield_curve.df(date2))
print(yield_curve.df(date3))
print(yield_curve.df(date4))
print(yield_curve.df(date5))

normalized_yield_curve = yield_curve.convert_to_convention("cont_act_365_curve", "Cont", "Act/365")
print("NORMALIZED YIELD CURVE:")
print(normalized_yield_curve.spot_rate(ref_date))
print(normalized_yield_curve.name)
print(normalized_yield_curve.base_date)
print(normalized_yield_curve.dates)
print(normalized_yield_curve.values)
print(normalized_yield_curve.curve)

date1 = Date.convert_str_to_datetime("2019-01-06")
date2 = Date.convert_str_to_datetime("2020-03-07")
print("Spot 2019-01-06: " + str(yield_curve.spot_rate(date1)))
print("Spot 2020-03-07: " + str(yield_curve.spot_rate(date2)))
print("Fwd rate: " + str(yield_curve.fwd_rate(date1, date2)))
print("Normalized spot 2019-01-06: " + str(normalized_yield_curve.spot_rate(date1)))
print("Normalized 2020-03-07: " + str(normalized_yield_curve.spot_rate(date2)))
print("Normalized fwd rate: " + str(normalized_yield_curve.fwd_rate(date1, date2)))
'''
