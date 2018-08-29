import datetime
import math
import Curve
from Data_Utils import DataUtils
from Date import Date
import collections
from Solver import Solver


class Bond(object):
    """Base class for bond objects"""

    _valid_coupon_frequency = ["ANNUAL", "A", "1", "SEMI-ANNUAL", "S", "2", "QUARTERLY", "Q", "4", "MONTHLY", "M", "12",
                               "DAILY", "D", "365"]
    _valid_yield_frequency = ["CONT", "CONTINUOUS", "ANNUAL", "A", "1", "SEMI-ANNUAL", "S", "2", "QUARTERLY", "Q", "4",
                              "MONTHLY", "M", "12", "DAILY", "D", "365", "SIMPLE", "SMP"]
    _valid_day_count_conv = ['1/1', 'ACT/ACT', 'ACT/365', 'ACT/360', '30/360', '30E/360', '30E/360 ISDA']

    _frequency_map = {
        "ANNUAL": {"inc": 1, "unit": "Y"},
        "A": {"inc": 1, "unit": "Y"},
        "1": {"inc": 1, "unit": "Y"},
        "SEMI-ANNUAL": {"inc": 6, "unit": "M"},
        "S": {"inc": 6, "unit": "M"},
        "2": {"inc": 6, "unit": "M"},
        "QUARTERLY": {"inc": 3, "unit": "M"},
        "Q": {"inc": 3, "unit": "M"},
        "4": {"inc": 3, "unit": "M"},
        "MONTHLY": {"inc": 1, "unit": "M"},
        "M": {"inc": 1, "unit": "M"},
        "12": {"inc": 1, "unit": "M"},
        "DAILY": {"inc": 1, "unit": "D"},
        "D": {"inc": 1, "unit": "D"},
        "365": {"inc": 1, "unit": "D"}
    }

    _compounding_periods = {
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

    def __init__(self, name, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                 issue_date=None, notional=100., clean_price=None, as_of_date=None, amort_schedule=None, oas=None,
                 discount_curve=None, ytm=None, yield_frequency=None, yield_day_count_conv=None, first_coupon_date=None,
                 last_coupon_date=None, reset_frequency=None, ref_rate_tenor=None, ref_rate_curve=None,
                 ref_rate_frequency=None, ref_rate_day_count_conv=None, initial_ref_rate=None, coupon_margin=None,
                 call_schedule=None, put_schedule=None, ytw=None, step_schedule=None, first_floating_date=None):

        self._name = name
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._coupon_frequency = coupon_frequency.upper()
        self._coupon_day_count_conv = coupon_day_count_conv.upper()
        self._dated_date = dated_date
        self._issue_date = issue_date
        self._notional = notional
        self._clean_price = clean_price
        self._as_of_date = as_of_date
        self._amort_schedule = amort_schedule
        self._oas = oas
        self._discount_curve = discount_curve
        self._ytm = ytm
        self._yield_frequency = yield_frequency
        self._yield_day_count_conv = yield_day_count_conv
        self._first_coupon_date = first_coupon_date
        self._last_coupon_date = last_coupon_date
        self._reset_frequency = reset_frequency
        self._ref_rate_tenor = ref_rate_tenor
        self._ref_rate_curve = ref_rate_curve
        self._ref_rate_frequency = ref_rate_frequency
        self._ref_rate_day_count_conv = ref_rate_day_count_conv
        self._initial_ref_rate = initial_ref_rate
        self._coupon_margin = coupon_margin
        self._call_schedule = call_schedule
        self._put_schedule = put_schedule
        self._ytw = ytw
        self._step_schedule = step_schedule
        self._first_floating_date = first_floating_date

        self.__validate_inputs()
        self.__complete_inputs()
        self._eff_first_coupon_date = Bond.__calc_eff_first_coupon_date(self._first_coupon_date, self._last_coupon_date,
                                                                        self._issue_date, self._dated_date,
                                                                        self._coupon_frequency)
        self._eff_last_coupon_date = Bond.__calc_eff_last_coupon_date(self._last_coupon_date, self._first_coupon_date,
                                                                      self._coupon_frequency, self._maturity_date)

        self._prin_cashflow_schedule = []
        self._prin_cashflow_values = []
        self._prin_schedule = self._generate_prin_schedule()
        self._prin_cashflows = self._generate_prin_cashflows()

        if self._reset_frequency is not None:

            self._reset_schedule = self._generate_int_cashflow_schedule(self._reset_frequency)
            self._reset_values = []

        self._resets = self._generate_resets()

        self._int_cashflow_schedule = self._generate_int_cashflow_schedule(self._coupon_frequency)
        self._int_cashflow_values = []
        self._int_cashflows = self._generate_int_cashflows()
        self._cashflow_schedule = []
        self._cashflow_values = []
        self._cashflows = self._generate_total_cashflows()

        self.__solve_yields_and_oas()
        self.__zero_out_analytics_post_maturity()

    @property
    def name(self):
        return self._name

    @property
    def maturity_date(self):
        return self._maturity_date

    @property
    def coupon(self):
        return self._coupon

    @property
    def dated_date(self):
        return self._dated_date

    @property
    def issue_date(self):
        return self._issue_date

    @property
    def oas(self):
        return self._oas

    @property
    def ytm(self):
        return self._ytm

    @property
    def ytw(self):
        return self._ytw

    @property
    def coupon_freq(self):
        return self._coupon_frequency

    @property
    def coupon_day_count(self):
        return self._coupon_day_count_conv

    @property
    def int_cashflow_schedule(self):
        return self._int_cashflow_schedule

    @property
    def int_cashflow_values(self):
        return self._int_cashflow_values

    @property
    def int_cashflows(self):
        return self._int_cashflows

    @property
    def prin_cashflow_schedule(self):
        return self._prin_cashflow_schedule

    @property
    def prin_cashflow_values(self):
        return self._prin_cashflow_values

    @property
    def prin_cashflows(self):
        return self._prin_cashflows

    @property
    def cashflow_schedule(self):
        return self._cashflow_schedule

    @property
    def cashflow_values(self):
        return self._cashflow_values

    @property
    def cashflows(self):
        return self._cashflows

    @staticmethod
    def __validate_inputs_lambda(name, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                                 issue_date, notional, clean_price, as_of_date, amort_schedule, oas, discount_curve,
                                 ytm, yield_frequency, yield_day_count_conv, first_coupon_date, last_coupon_date,
                                 reset_frequency, ref_rate_tenor, ref_rate_curve, ref_rate_frequency,
                                 ref_rate_day_count_conv, initial_ref_rate, coupon_margin, call_schedule, put_schedule,
                                 ytw, step_schedule, first_floating_date):

        DataUtils.validate_type(name, str, "Name")

        if coupon is not None:
            DataUtils.validate_number(coupon, "Coupon rate")

        DataUtils.validate_type(maturity_date, datetime.date, "Maturity date")
        DataUtils.validate_in_list(coupon_frequency, Bond._valid_coupon_frequency, "Coupon frequency")
        DataUtils.validate_in_list(coupon_day_count_conv, Bond._valid_day_count_conv, "Coupon day count convention")
        Bond.__validate_start_date(dated_date, "Dated date", maturity_date)

        if issue_date is not None:
            Bond.__validate_start_date(issue_date, "Issue date", maturity_date)
        else:
            issue_date = dated_date

        DataUtils.validate_number(notional, "Notional")

        if clean_price is not None:
            DataUtils.validate_number(clean_price, "Clean price")

        if as_of_date is not None:
            DataUtils.validate_type(as_of_date, datetime.date, "As of date")

        if amort_schedule is not None:
            Bond.__validate_amort_or_step_schedule(amort_schedule, "Amortization", issue_date, maturity_date)

        if oas is not None:
            DataUtils.validate_number(oas, "OAS")

        if discount_curve is not None:
            DataUtils.validate_type(discount_curve, Curve.YieldCurve, "Discount curve")

        if ytm is not None:
            DataUtils.validate_number(ytm, "YTM")

        if yield_frequency is not None:

            yield_frequency = yield_frequency.upper()
            DataUtils.validate_in_list(yield_frequency, Bond._valid_yield_frequency, "Yield compounding frequency")

        if yield_day_count_conv is not None:

            yield_day_count_conv = yield_day_count_conv.upper()
            DataUtils.validate_in_list(yield_day_count_conv, Bond._valid_day_count_conv, "Yield day count convention")

        if first_coupon_date is not None:
            Bond.__validate_coupon_date(first_coupon_date, "First coupon date", maturity_date, issue_date, dated_date,
                                        first_coupon_date, last_coupon_date)

        if last_coupon_date is not None:
            Bond.__validate_coupon_date(last_coupon_date, "Last coupon date", maturity_date, issue_date, dated_date,
                                        first_coupon_date, last_coupon_date)

        if reset_frequency is not None:
            DataUtils.validate_in_list(reset_frequency, Bond._valid_coupon_frequency, "Reset frequency")

        if ref_rate_tenor is not None:
            DataUtils.validate_tenor(ref_rate_tenor)

        if ref_rate_curve is not None:
            DataUtils.validate_type(ref_rate_curve, Curve.YieldCurve, "Reference rate curve")

        if ref_rate_frequency is not None:
            DataUtils.validate_in_list(ref_rate_frequency, Bond._valid_coupon_frequency, "Reference rate frequency")

        if ref_rate_day_count_conv is not None:
            DataUtils.validate_in_list(ref_rate_day_count_conv, Bond._valid_day_count_conv,
                                       "Reference rate day count convention")

        if initial_ref_rate is not None:
            DataUtils.validate_number(initial_ref_rate, "Initial reference rate")

        if coupon_margin is not None:
            DataUtils.validate_number(initial_ref_rate, "Coupon margin")

        if call_schedule is not None:
            Bond.__validate_call_put_schedule(call_schedule, issue_date, maturity_date)

        if put_schedule is not None:
            Bond.__validate_call_put_schedule(put_schedule, issue_date, maturity_date)

        if ytw is not None:
            DataUtils.validate_number(ytm, "YTW")

        if step_schedule is not None:
            Bond.__validate_amort_or_step_schedule(step_schedule, "Step", issue_date, maturity_date)

        if first_floating_date is not None:
            Bond.__validate_coupon_date(first_floating_date, "First reset date", maturity_date, issue_date, dated_date,
                                        first_coupon_date, last_coupon_date)

    @staticmethod
    def __quick_validate_inputs(name, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                                issue_date, notional, clean_price, as_of_date, amort_schedule, oas, discount_curve, ytm,
                                yield_frequency, yield_day_count_conv, first_coupon_date, last_coupon_date,
                                reset_frequency, ref_rate_tenor, ref_rate_curve, ref_rate_frequency,
                                ref_rate_day_count_conv, initial_ref_rate, coupon_margin, call_schedule, put_schedule,
                                ytw, step_schedule, first_floating_date):

        Bond.__validate_inputs_lambda(name, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                                      issue_date, notional, clean_price, as_of_date, amort_schedule, oas,
                                      discount_curve, ytm, yield_frequency, yield_day_count_conv, first_coupon_date,
                                      last_coupon_date, reset_frequency, ref_rate_tenor, ref_rate_curve,
                                      ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate, coupon_margin,
                                      call_schedule, put_schedule, ytw, step_schedule, first_floating_date)

    def __validate_inputs(self):
        Bond.__validate_inputs_lambda(self._name, self._maturity_date, self._coupon, self._coupon_frequency,
                                      self._coupon_day_count_conv, self._dated_date, self._issue_date,
                                      self._notional, self._clean_price, self._as_of_date, self._amort_schedule,
                                      self._oas, self._discount_curve, self._ytm, self._yield_frequency,
                                      self._yield_day_count_conv, self._first_coupon_date, self._last_coupon_date,
                                      self._reset_frequency, self._ref_rate_tenor, self._ref_rate_curve,
                                      self._ref_rate_frequency, self._ref_rate_day_count_conv, self._initial_ref_rate,
                                      self._coupon_margin, self._call_schedule, self._put_schedule, self._ytw,
                                      self._step_schedule, self._first_floating_date)

    def __complete_inputs(self):

        if self._issue_date is None:
            self._issue_date = self._dated_date

        self._notional = float(self._notional)

        if self._clean_price is not None:
            self._clean_price = float(self._clean_price)

        if self._oas is not None:
            self._oas = float(self._oas)

        if self._discount_curve is not None:
            self._discount_curve = self._discount_curve.convert_to_convention("normalized_curve",
                                                                              self._coupon_frequency,
                                                                              self._coupon_day_count_conv)

        if self._ytm is not None:
            self._ytm = float(self._ytm)

        if self._yield_frequency is not None:
            self._yield_frequency = self._yield_frequency.upper()
        else:
            self._yield_frequency = self._coupon_frequency

        if self._yield_day_count_conv is not None:
            self._yield_day_count_conv = self._yield_day_count_conv.upper()
        else:
            self._yield_day_count_conv = self._coupon_day_count_conv

        if (self._ref_rate_frequency is None) and (self._ref_rate_curve is not None):
            self._ref_rate_frequency = self._reset_frequency

        if (self._ref_rate_day_count_conv is None) and (self._ref_rate_curve is not None):
            self._ref_rate_day_count_conv = self._coupon_day_count_conv

        if self._coupon_margin is not None:
            self._coupon_margin = float(self._coupon_margin)

        if self._ytw is not None:
            self._ytw = float(self._ytw)

    @staticmethod
    def __validate_start_date(date, name, maturity_date):

        DataUtils.validate_type(date, datetime.date, name)

        if date >= maturity_date:
            raise ValueError(name + " must be less than maturity date")

    @staticmethod
    def __validate_coupon_date(date, name, maturity_date, issue_date, dated_date, first_coupon_date, last_coupon_date):

        DataUtils.validate_type(date, datetime.date, name)

        if (date > maturity_date) or (date < issue_date) or (date < dated_date):
            raise ValueError(name + " must be less than or equal to maturity date, "
                             "and greater than or equal to both dated date and issue date")

        if (first_coupon_date is not None) and (last_coupon_date is not None) and \
                (first_coupon_date >= last_coupon_date):
            raise ValueError("First coupon date must be less than last coupon date")

    @staticmethod
    def __validate_pricing_input(date, ytm, oas, discount_curve, ref_rate_curve):

        if (ytm is None) and ((oas is None) or (discount_curve is None)):
            raise ValueError("Bond pricing requires either YTM, or OAS and discount curve")

        if (ytm is None) and (discount_curve.base_date != date):
            raise ValueError("Discount curve base date must equal valuation date")

        if (ref_rate_curve is not None) and (ref_rate_curve.base_date != date):
            raise ValueError("Reference rate curve base date must equal valuation date")

    @staticmethod
    def __validate_amort_or_step_schedule(amort_schedule, name, issue_date, maturity_date):

        DataUtils.validate_type(amort_schedule, collections.OrderedDict, name + " schedule")

        dates = list(amort_schedule.keys())
        values = list(amort_schedule.values())

        if len(dates) == 0:
            raise ValueError(name + " schedule must be of length 1 or greater")

        [DataUtils.validate_type(date, datetime.date, name + " date") for date in dates]
        [DataUtils.validate_number(val, name + " value") for val in values]

        if dates[0] < issue_date:
            raise ValueError("First date must be greater than or equal to issue date")

        if dates[len(dates) - 1] > maturity_date:
            raise ValueError("Last date must be less than or equal to maturity date")

        if len(dates) == 1:
            return

        if not all(x < y for x, y in zip(dates, dates[1:])):
            raise ValueError(name + " dates must be strictly increasing")

    @staticmethod
    def __validate_call_put_schedule(schedule, issue_date, maturity_date):

        DataUtils.validate_type(schedule, collections.OrderedDict, "Call or put schedule")

        for date, val in schedule.items():

            DataUtils.validate_type(date, datetime.date, "Date")
            DataUtils.validate_type(val, dict, "Strike and exercise type")
            DataUtils.validate_number(val["Strike"], "Strike")
            DataUtils.validate_in_list(val["Exercise Type"], ["American", "European", "Bermudan"], "Exercise type")

            if (date < issue_date) or (date >= maturity_date):
                raise ValueError("Call or put date must be greater than or equal to issue date "
                                 "and less than maturity date")

    @staticmethod
    def __calc_eff_first_coupon_date(first_coupon_date, last_coupon_date, issue_date, dated_date, coupon_frequency):

        if first_coupon_date:
            return first_coupon_date

        elif last_coupon_date:

            base_date = last_coupon_date
            iter_date = base_date
            later_of_issue_or_dated_date = issue_date if (issue_date >= dated_date) else dated_date
            unit_count = 0
            freq_dict = Bond._frequency_map[coupon_frequency]

            while iter_date >= later_of_issue_or_dated_date:

                unit_count += freq_dict["inc"]
                iter_date = Date.add_tenor(base_date, "-" + str(unit_count) + str(freq_dict["unit"]))

            iter_date = Date.add_tenor(iter_date, str(freq_dict["inc"]) + str(freq_dict["unit"]))
            return iter_date

        else:
            return issue_date if (issue_date >= dated_date) else dated_date

    @staticmethod
    def __calc_eff_last_coupon_date(last_coupon_date, first_coupon_date, coupon_frequency, maturity_date):

        if last_coupon_date:
            return last_coupon_date

        elif first_coupon_date:

            base_date = first_coupon_date
            iter_date = base_date
            unit_count = 0
            freq_dict = Bond._frequency_map[coupon_frequency]

            while iter_date <= maturity_date:

                unit_count += freq_dict["inc"]
                iter_date = Date.add_tenor(base_date, str(unit_count) + str(freq_dict["unit"]))

            iter_date = Date.add_tenor(iter_date, "-" + str(freq_dict["inc"]) + str(freq_dict["unit"]))
            return iter_date

        else:
            return maturity_date

    @staticmethod
    def __generate_int_cashflow_schedule_lambda(dated_date, coupon_frequency, first_coupon_date, eff_first_coupon_date,
                                                eff_last_coupon_date):

        cf_schedule = [dated_date]

        if first_coupon_date:

            base_date = eff_first_coupon_date
            iter_date = base_date
            unit_count = 0
            freq_dict = Bond._frequency_map[coupon_frequency]

            while iter_date < eff_last_coupon_date:

                if iter_date > dated_date:
                    cf_schedule.append(iter_date)

                unit_count += freq_dict["inc"]
                iter_date = Date.add_tenor(base_date, str(unit_count) + str(freq_dict["unit"]))

            if eff_last_coupon_date > dated_date:
                cf_schedule.append(eff_last_coupon_date)

        else:

            base_date = eff_last_coupon_date
            iter_date = base_date
            unit_count = 0
            freq_dict = Bond._frequency_map[coupon_frequency]

            while iter_date > eff_first_coupon_date:
                cf_schedule.insert(1, iter_date)

                unit_count += freq_dict["inc"]
                iter_date = Date.add_tenor(base_date, "-" + str(unit_count) + str(freq_dict["unit"]))

            if eff_first_coupon_date > dated_date:
                cf_schedule.insert(1, eff_first_coupon_date)

        return cf_schedule

    @staticmethod
    def __quick_generate_int_cashflow_schedule(dated_date, coupon_frequency, first_coupon_date, eff_first_coupon_date,
                                               eff_last_coupon_date):

        return Bond.__generate_int_cashflow_schedule_lambda(dated_date, coupon_frequency, first_coupon_date,
                                                            eff_first_coupon_date, eff_last_coupon_date)

    def _generate_int_cashflow_schedule(self, frequency):

        if self._first_floating_date is None:

            return Bond.__generate_int_cashflow_schedule_lambda(self._dated_date, frequency, self._first_coupon_date,
                                                                self._eff_first_coupon_date, self._eff_last_coupon_date)

        else:

            return Bond.__generate_int_cashflow_schedule_lambda(self._dated_date, frequency, self._first_floating_date,
                                                                self._first_floating_date, self._eff_last_coupon_date)

    @staticmethod
    def _generate_resets_lambda(reset_schedule, ref_rate_curve, ref_rate_tenor, initial_ref_rate, first_floating_date,
                                ref_rate_frequency, ref_rate_day_count_conv):

        resets = collections.OrderedDict()
        base_date = ref_rate_curve.base_date
        normalized_ref_rate_curve = ref_rate_curve.convert_to_convention("normalized_ref_rate_curve",
                                                                         ref_rate_frequency, ref_rate_day_count_conv)

        for i in range(0, len(reset_schedule)):

            date = reset_schedule[i]

            if i == 0:
                prev_date = date
            else:
                prev_date = reset_schedule[i - 1]

            if date < base_date:
                reset_rate = 0.
            else:

                reset_rate = normalized_ref_rate_curve.fwd_rate(date, Date.add_tenor(date, ref_rate_tenor))

                if (first_floating_date is None) and ((i == 0) or (prev_date < base_date)):
                    resets[prev_date] = initial_ref_rate

            resets[date] = reset_rate

        return resets

    @staticmethod
    def _generate_int_cashflows_lambda(int_schedule, prin_schedule, coupon, resets, coupon_margin,
                                       eff_first_coupon_date, eff_last_coupon_date, coupon_day_count_conv, calc_coupon,
                                       step_schedule, first_floating_date):

        int_cashflows = collections.OrderedDict()
        int_cashflows[int_schedule[0]] = 0.

        for i in range(1, len(int_schedule)):

            date = int_schedule[i]

            if (date >= eff_first_coupon_date) and (date <= eff_last_coupon_date):

                outstanding_principal = Bond._get_schedule_val(prin_schedule, date)
                year_frac = Date.year_frac(int_schedule[i - 1], date, coupon_day_count_conv)

                if (coupon is None) and (resets is None):
                    cashflow_val = calc_coupon(date) * year_frac
                elif resets is None:
                    cashflow_val = calc_coupon(date, outstanding_principal, coupon, step_schedule) * year_frac
                elif first_floating_date is None:
                    cashflow_val = calc_coupon(date, outstanding_principal, resets, coupon_margin) * year_frac
                else:
                    cashflow_val = calc_coupon(date, outstanding_principal, resets, coupon_margin,
                                               first_floating_date) * year_frac

                int_cashflows[date] = cashflow_val

        return int_cashflows

    @staticmethod
    def _generate_prin_cashflows_lambda(notional, maturity_date, amort_schedule, is_payment):

        prin_cashflows = collections.OrderedDict()

        if amort_schedule is None:
            prin_cashflows[maturity_date] = notional
        else:

            outstanding_principal = notional
            dates = list(amort_schedule.keys())

            for i in range(0, len(dates)):

                date = dates[i]

                prin_payment = outstanding_principal * (amort_schedule[date] * 0.01)
                prin_cashflows[date] = prin_payment if is_payment else outstanding_principal
                outstanding_principal -= prin_payment

        return prin_cashflows

    @staticmethod
    def _generate_total_cashflows_lambda(prin_cashflows, int_cashflows):

        total_cashflows = collections.OrderedDict(int_cashflows)
        int_cashflow_dates = int_cashflows.keys()

        for date in prin_cashflows.keys():

            if date in int_cashflow_dates:
                total_cashflows[date] += prin_cashflows[date]
            else:
                total_cashflows[date] = prin_cashflows[date]

        total_cashflows = collections.OrderedDict(sorted(total_cashflows.items(), key=lambda x: x[0]))

        return total_cashflows

    @classmethod
    def _quick_generate_resets(cls, reset_schedule, ref_rate_curve, ref_rate_tenor, initial_ref_rate,
                               first_floating_date, ref_rate_frequency, ref_rate_day_count_conv):

        return Bond._generate_resets_lambda(reset_schedule, ref_rate_curve, ref_rate_tenor, initial_ref_rate,
                                            first_floating_date, ref_rate_frequency, ref_rate_day_count_conv)

    @classmethod
    def _quick_generate_int_cashflows(cls, int_schedule, prin_schedule, coupon, resets, coupon_margin,
                                      eff_first_coupon_date, eff_last_coupon_date, coupon_day_count_conv,
                                      step_schedule, first_floating_date):

        return Bond._generate_int_cashflows_lambda(int_schedule, prin_schedule, coupon, resets, coupon_margin,
                                                   eff_first_coupon_date, eff_last_coupon_date, coupon_day_count_conv,
                                                   cls._calc_coupon_lambda, step_schedule, first_floating_date)

    @classmethod
    def _quick_generate_prin_cashflows(cls, notional, maturity_date, amort_schedule):
        return Bond._generate_prin_cashflows_lambda(notional, maturity_date, amort_schedule, True)

    @classmethod
    def _quick_generate_prin_schedule(cls, notional, maturity_date, amort_schedule):
        return Bond._generate_prin_cashflows_lambda(notional, maturity_date, amort_schedule, False)

    @classmethod
    def _quick_generate_total_cashflows(cls, prin_cashflows, int_cashflows):
        return Bond._generate_total_cashflows_lambda(prin_cashflows, int_cashflows)

    def _generate_resets(self):

        if self._ref_rate_curve is None:
            return None

        resets = Bond._generate_resets_lambda(self._reset_schedule, self._ref_rate_curve, self._ref_rate_tenor,
                                              self._initial_ref_rate, self._first_floating_date,
                                              self._ref_rate_frequency, self._ref_rate_day_count_conv)

        self._reset_values = resets.values()

        return resets

    def _generate_int_cashflows(self):

        int_cashflows = Bond._generate_int_cashflows_lambda(self._int_cashflow_schedule, self._prin_schedule, None,
                                                            None, None, self._eff_first_coupon_date,
                                                            self._eff_last_coupon_date, self._coupon_day_count_conv,
                                                            Bond._calc_coupon_lambda, self._step_schedule,
                                                            self._first_floating_date)

        self._int_cashflow_values = int_cashflows.values()

        return int_cashflows

    def _generate_prin_schedule(self):

        prin_schedule = Bond._generate_prin_cashflows_lambda(self._notional, self._maturity_date, self._amort_schedule,
                                                             False)

        return prin_schedule

    def _generate_prin_cashflows(self):

        prin_cashflows = Bond._generate_prin_cashflows_lambda(self._notional, self._maturity_date, self._amort_schedule,
                                                              True)

        self._prin_cashflow_schedule = prin_cashflows.keys()
        self._prin_cashflow_values = prin_cashflows.values()

        return prin_cashflows

    def _generate_total_cashflows(self):

        total_cashflows = Bond._generate_total_cashflows_lambda(self._prin_cashflows, self._int_cashflows)

        self._total_cashflow_schedule = total_cashflows.keys()
        self._total_cashflow_values = total_cashflows.values()

        return total_cashflows

    @staticmethod
    def _calc_coupon_lambda(date):
        return 0.

    def dirty_price(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        Bond.__validate_pricing_input(val_date, self._ytm, self._oas, self._discount_curve, self._ref_rate_curve)

        price = 0.

        if (self._oas is not None) and (self._discount_curve is not None):

            discount_curve_with_oas = self._discount_curve.parallel_shift(self._oas)
            df = discount_curve_with_oas.df
            fwd_df = discount_curve_with_oas.fwd_df

        else:

            df = self.__ytm_df(val_date, self._ytm, self._yield_day_count_conv, self._yield_frequency,
                               self._maturity_date)
            fwd_df = self.__ytm_fwd_df(self._ytm, self._yield_day_count_conv, self._yield_frequency,
                                       self._maturity_date)

        oa_cashflows = Bond.__calc_oa_cashflows(val_date, self._cashflows, fwd_df, self._call_schedule,
                                                self._put_schedule)

        for cf_date, cf_val in oa_cashflows.items():

            if cf_date > val_date:
                price += cf_val * df(cf_date)

        return price

    def __zero_out_analytics_post_maturity(self):

        if (self._as_of_date is not None) and (self._as_of_date >= self._maturity_date):

            self._ytm = 0.
            self._oas = 0.
            self._clean_price = 0.

    @staticmethod
    def __ytm_df(date1, ytm, yield_day_count_conv, yield_frequency, maturity_date=None):

        def df(date2):
            return Bond.__calc_df(date1, date2, ytm, yield_day_count_conv, yield_frequency, maturity_date)

        return df

    @staticmethod
    def __ytm_fwd_df(ytm, yield_day_count_conv, yield_frequency, maturity_date=None):

        def fwd_df(date1, date2):
            return Bond.__calc_df(date1, date2, ytm, yield_day_count_conv, yield_frequency, maturity_date)

        return fwd_df

    @staticmethod
    def __calc_df(date1, date2, ytm, yield_day_count_conv, yield_frequency, maturity_date=None):

        DataUtils.validate_type(date1, datetime.date, "Date")
        DataUtils.validate_type(date2, datetime.date, "Date")

        if date1 > date2:
            raise ValueError("End date must be greater than or equal to start date")

        if yield_day_count_conv != "30E/360 ISDA":
            year_frac = Date.year_frac(date1, date2, yield_day_count_conv)
        else:
            year_frac = Date.year_frac(date1, date2, yield_day_count_conv, maturity_date)

        if (yield_frequency == "CONT") or (yield_frequency == "CONTINUOUS"):
            return math.exp(-ytm * year_frac)

        elif (yield_frequency == "SMP") or (yield_frequency == "SIMPLE"):
            return 1 / (1 + ytm * year_frac)

        else:

            n = Bond._compounding_periods[yield_frequency]

            if (1 + (ytm / n)) == 0.:
                return 0.

            return (1 + (ytm / n)) ** (-year_frac * n)

    @staticmethod
    def _get_schedule_val(schedule, val):

        dates = list(schedule.keys())

        if len(dates) == 0:
            return None

        if val < dates[0]:
            return schedule[dates[0]]

        if val > dates[len(dates) - 1]:
            return schedule[dates[len(dates) - 1]]

        if val in dates:
            return schedule[val]

        for i in range(1, len(dates)):

            if val < dates[i]:
                return schedule[dates[i]]

    @staticmethod
    def __acc_int_lambda(val_date, dated_date, eff_last_coupon_date, int_cashflows, coupon_day_count_conv):

        DataUtils.validate_type(val_date, datetime.date, "Valuation date")

        cf_dates = list(int_cashflows.keys())
        period_start_date = None
        period_end_date = None
        cf_val = None

        if (val_date <= dated_date) or (val_date >= eff_last_coupon_date) or (val_date >= cf_dates[len(cf_dates) - 1]):
            return 0.

        for cf_date in cf_dates:

            if val_date < cf_date:

                period_start_date = cf_dates[cf_dates.index(cf_date) - 1]
                period_end_date = cf_date

                cf_val = int_cashflows[period_end_date]

                break

        return cf_val * Date.partial_frac(period_start_date, val_date, period_end_date, coupon_day_count_conv)

    @staticmethod
    def __quick_acc_int(val_date, dated_date, eff_last_coupon_date, int_cashflows, coupon_day_count_conv):
        return Bond.__acc_int_lambda(val_date, dated_date, eff_last_coupon_date, int_cashflows, coupon_day_count_conv)

    def acc_int(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        return Bond.__acc_int_lambda(val_date, self._dated_date, self._eff_last_coupon_date, self._int_cashflows,
                                     self._coupon_day_count_conv)

    def clean_price(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        if (self._clean_price is not None) and (val_date == self._as_of_date):
            return self._clean_price

        dirty_price = self.dirty_price(val_date)
        ai = self.acc_int(val_date)

        return dirty_price - ai

    @staticmethod
    def __price_lambda(val_date, dirty_or_clean, maturity_date, coupon_frequency, coupon_day_count_conv,
                       dated_date, cashflows, int_cashflows, issue_date=None, notional=100., amort_schedule=None,
                       oas=None, discount_curve=None, ytm=None, yield_frequency=None, yield_day_count_conv=None,
                       first_coupon_date=None, last_coupon_date=None, reset_frequency=None, ref_rate_tenor=None,
                       ref_rate_curve=None, ref_rate_frequency=None, ref_rate_day_count_conv=None,
                       initial_ref_rate=None, coupon_margin=None, call_schedule=None, put_schedule=None, ytw=None,
                       step_schedule=None, first_floating_date=None):

        if (val_date >= maturity_date) or (len(cashflows.keys()) == 0):
            return 0.

        Bond.__quick_validate_inputs("", maturity_date, None, coupon_frequency, coupon_day_count_conv, dated_date,
                                     issue_date, notional, 0., val_date, amort_schedule, oas, discount_curve, ytm,
                                     yield_frequency, yield_day_count_conv, first_coupon_date, last_coupon_date,
                                     reset_frequency, ref_rate_tenor, ref_rate_curve, ref_rate_frequency,
                                     ref_rate_day_count_conv, initial_ref_rate, coupon_margin, call_schedule,
                                     put_schedule, ytw, step_schedule, first_floating_date)

        dirty_or_clean = dirty_or_clean.upper()
        DataUtils.validate_in_list(dirty_or_clean, ["DIRTY", "CLEAN"], "'Dirty or clean'")

        if oas is not None:
            oas = float(oas)

        if discount_curve is not None:
            discount_curve = discount_curve.convert_to_convention("normalized_curve", coupon_frequency,
                                                                  coupon_day_count_conv)

        if ytm is not None:
            ytm = float(ytm)

        if yield_frequency is not None:
            yield_frequency = yield_frequency.upper()
        else:
            yield_frequency = coupon_frequency

        if yield_day_count_conv is not None:
            yield_day_count_conv = yield_day_count_conv.upper()
        else:
            yield_day_count_conv = coupon_day_count_conv

        Bond.__validate_pricing_input(val_date, ytm, oas, discount_curve, ref_rate_curve)

        eff_last_coupon_date = Bond.__calc_eff_last_coupon_date(last_coupon_date, first_coupon_date, coupon_frequency,
                                                                maturity_date)

        if (oas is not None) and (discount_curve is not None):

            discount_curve_with_oas = discount_curve.parallel_shift(oas)
            df = discount_curve_with_oas.df
            fwd_df = discount_curve_with_oas.fwd_df

        else:

            df = Bond.__ytm_df(val_date, ytm, yield_day_count_conv, yield_frequency, maturity_date)
            fwd_df = Bond.__ytm_fwd_df(ytm, yield_day_count_conv, yield_frequency, maturity_date)

        price = Bond.__calc_price(val_date, cashflows, df, fwd_df, call_schedule, put_schedule)

        if dirty_or_clean == "CLEAN":

            ai = Bond.__quick_acc_int(val_date, dated_date, eff_last_coupon_date, int_cashflows, coupon_day_count_conv)
            price -= ai

        return price

    @staticmethod
    def __calc_oa_cashflows(val_date, cashflows, fwd_df, call_schedule, put_schedule):

        oa_cashflows = collections.OrderedDict(cashflows)

        if (call_schedule is not None) or (put_schedule is not None):

            def update_cashflows(prime_date, iter_dt, oa_cfs):

                breached_val_dt = False
                future_cashflows = collections.OrderedDict({k: v for k, v in oa_cfs.items() if k > iter_dt})
                pv = 0.

                for cf_date, cf_val in future_cashflows.items():

                    if cf_date <= val_date:
                        breached_val_dt = True
                        break

                    if cf_date > iter_dt:
                        pv += cf_val * fwd_df(iter_dt, cf_date)

                call_intrinsic_value = 0
                put_intrinsic_value = 0
                call_strike = 0.
                put_strike = 0.

                if prime_date in call_dates:

                    call_strike = call_schedule[prime_date]["Strike"]

                    if pv > call_strike:
                        call_intrinsic_value = pv - call_strike

                if prime_date in put_dates:

                    put_strike = put_schedule[prime_date]["Strike"]

                    if pv < put_strike:
                        put_intrinsic_value = put_strike - pv

                if (call_intrinsic_value == 0) and (put_intrinsic_value == 0):
                    return [oa_cfs, breached_val_dt]
                else:
                    strike = 0.

                    if call_intrinsic_value > 0:
                        strike = call_strike
                    elif put_intrinsic_value > 0:
                        strike = put_strike

                    if iter_dt not in oa_cfs:
                        oa_cfs[iter_dt] = strike
                    else:
                        oa_cfs[iter_dt] += strike

                oa_cfs = collections.OrderedDict(sorted(oa_cfs.items(), key=lambda x: x[0]))
                cashflow_dates = list(oa_cfs.keys())

                for j in range(cashflow_dates.index(iter_dt) + 1, len(cashflow_dates)):
                    del oa_cfs[cashflow_dates[j]]

                return [oa_cfs, breached_val_dt]

            tmp_cf_dates = list(oa_cashflows.keys())
            call_dates = list(call_schedule.keys()) if (call_schedule is not None) else []
            put_dates = list(put_schedule.keys()) if (put_schedule is not None) else []

            if (call_schedule is not None) and (put_schedule is not None):
                exercise_dates = sorted(list(set(call_dates + put_dates)))
            elif call_schedule is not None:
                exercise_dates = call_dates
            else:
                exercise_dates = put_dates

            breached_val_date = False

            for date in reversed(exercise_dates):

                call_option_type = ""
                put_option_type = ""

                if date in call_dates:
                    call_option_type = call_schedule[date]["Exercise Type"]

                if date in put_dates:
                    put_option_type = put_schedule[date]["Exercise Type"]

                american_in_effect = (call_option_type == "American") or (put_option_type == "American")

                if (date >= val_date) or american_in_effect:

                    if breached_val_date:
                        break

                    if american_in_effect:

                        if (call_option_type == "American") and (put_option_type == "American"):

                            call_date_idx = call_dates.index(date)
                            put_date_idx = put_dates.index(date)

                            if (call_date_idx == (len(call_dates) - 1)) or (put_date_idx == (len(put_dates) - 1)):
                                iter_date = tmp_cf_dates[-1] - datetime.timedelta(1)
                            else:

                                call_date = call_dates[call_date_idx]
                                next_call_date = call_dates[call_date_idx + 1]
                                call_duration = next_call_date - call_date

                                put_date = call_dates[put_date_idx]
                                next_put_date = put_dates[put_date_idx + 1]
                                put_duration = next_put_date - put_date

                                if call_duration >= put_duration:
                                    iter_date = next_call_date - datetime.timedelta(1)
                                else:
                                    iter_date = next_put_date - datetime.timedelta(1)

                        elif call_option_type == "American":

                            date_idx = call_dates.index(date)

                            if date_idx == (len(call_dates) - 1):
                                iter_date = tmp_cf_dates[-1] - datetime.timedelta(1)
                            else:
                                iter_date = call_dates[date_idx + 1] - datetime.timedelta(1)

                        else:

                            date_idx = put_dates.index(date)

                            if date_idx == (len(put_dates) - 1):
                                iter_date = tmp_cf_dates[-1] - datetime.timedelta(1)
                            else:
                                iter_date = put_dates[date_idx + 1] - datetime.timedelta(1)

                        while iter_date >= date:

                            if breached_val_date:
                                break

                            results = update_cashflows(date, iter_date, oa_cashflows)
                            oa_cashflows = results[0]
                            breached_val_date = results[1]

                            iter_date -= datetime.timedelta(1)

                    else:

                        results = update_cashflows(date, date, oa_cashflows)
                        oa_cashflows = results[0]
                        breached_val_date = results[1]

                else:
                    break

        return oa_cashflows

    @staticmethod
    def __calc_price(val_date, cashflows, df, fwd_df, call_schedule, put_schedule):

        oa_cashflows = Bond.__calc_oa_cashflows(val_date, cashflows, fwd_df, call_schedule, put_schedule)
        price = 0.

        for cf_date, cf_val in oa_cashflows.items():

            if cf_date > val_date:
                price += cf_val * df(cf_date)

        return price

    @classmethod
    def quick_price(cls, val_date, dirty_or_clean, maturity_date, coupon, coupon_frequency, coupon_day_count_conv,
                    dated_date, issue_date=None, notional=100., amort_schedule=None, oas=None, discount_curve=None,
                    ytm=None, yield_frequency=None, yield_day_count_conv=None, first_coupon_date=None,
                    last_coupon_date=None, reset_frequency=None, ref_rate_tenor=None, ref_rate_curve=None,
                    ref_rate_frequency=None, ref_rate_day_count_conv=None, initial_ref_rate=None, coupon_margin=None,
                    call_schedule=None, put_schedule=None, ytw=None, step_schedule=None, first_floating_date=None):

        if val_date >= maturity_date:
            return 0.

        Bond.__quick_validate_inputs("", maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                                     issue_date, notional, 0., val_date, amort_schedule, oas, discount_curve, ytm,
                                     yield_frequency, yield_day_count_conv, first_coupon_date, last_coupon_date,
                                     reset_frequency, ref_rate_tenor, ref_rate_curve, ref_rate_frequency,
                                     ref_rate_day_count_conv, initial_ref_rate, coupon_margin, call_schedule,
                                     put_schedule, ytw, step_schedule, first_floating_date)

        if issue_date is None:
            issue_date = dated_date

        notional = float(notional)

        eff_first_coupon_date = Bond.__calc_eff_first_coupon_date(first_coupon_date, last_coupon_date, issue_date,
                                                                  dated_date, coupon_frequency)
        eff_last_coupon_date = Bond.__calc_eff_last_coupon_date(last_coupon_date, first_coupon_date, coupon_frequency,
                                                                maturity_date)

        prin_cashflows = cls._quick_generate_prin_cashflows(notional, maturity_date, amort_schedule)
        prin_schedule = cls._quick_generate_prin_schedule(notional, maturity_date, amort_schedule)

        if ref_rate_curve is not None:

            reset_schedule = Bond.__quick_generate_int_cashflow_schedule(dated_date, reset_frequency, first_coupon_date,
                                                                         eff_first_coupon_date, eff_last_coupon_date)
            resets = cls._quick_generate_resets(reset_schedule, ref_rate_curve, ref_rate_tenor, initial_ref_rate,
                                                first_floating_date, ref_rate_frequency, ref_rate_day_count_conv)

        else:
            resets = None

        int_cashflow_schedule = Bond.__quick_generate_int_cashflow_schedule(dated_date, coupon_frequency,
                                                                            first_coupon_date, eff_first_coupon_date,
                                                                            eff_last_coupon_date)
        int_cashflows = cls._quick_generate_int_cashflows(int_cashflow_schedule, prin_schedule, coupon, resets,
                                                          coupon_margin, eff_first_coupon_date, eff_last_coupon_date,
                                                          coupon_day_count_conv, step_schedule, first_floating_date)
        total_cashflows = cls._quick_generate_total_cashflows(prin_cashflows, int_cashflows)

        return cls.__price_lambda(val_date, dirty_or_clean, maturity_date, coupon_frequency, coupon_day_count_conv,
                                  dated_date, total_cashflows, int_cashflows, issue_date, notional, oas, discount_curve,
                                  ytm, yield_frequency, yield_day_count_conv, first_coupon_date, last_coupon_date,
                                  call_schedule, put_schedule, ytw, step_schedule)

    def __solve_yields_and_oas(self):

        if (self._ytm is None) and (self._as_of_date is not None) and ((self._clean_price is not None) or
                                                                       ((self._oas is not None) and
                                                                        (self._discount_curve is not None))):

            self._ytm = self.__solve_yield("YTM")

        if (self._oas is None) and (self._discount_curve is not None) and (self._as_of_date is not None) and \
                ((self._clean_price is not None) or (self._ytm is not None)):

            self._oas = self.__solve_yield("OAS")

        if (self._ytw is None) and (self._as_of_date is not None) and ((self._call_schedule is not None) or
                                                                       (self._put_schedule is not None)) and \
                ((self._clean_price is not None) or ((self._oas is not None) and (self._discount_curve is not None))):

            self._ytw = self.__solve_yield("YTW")

        elif (self._call_schedule is None) and (self._put_schedule is None):
            self._ytw = self._ytm

    def __solve_yield(self, ytm_oas_or_ytw):

        if self._as_of_date >= self._maturity_date:
            return 0.

        ytm_oas_or_ytw = ytm_oas_or_ytw.upper()
        DataUtils.validate_in_list(ytm_oas_or_ytw, ["YTM", "OAS", "YTW"], "'YTM, OAS, or YTW'")

        a = -1.
        b = 1.
        max_iter = 1000000
        target_tol = 0.0000000001
        interval_tol = 0.00000000000000000001

        iter_price = self.__get_iter_price(ytm_oas_or_ytw)

        return Solver.bisection(iter_price, a, b, target_tol, interval_tol, max_iter)

    def __get_iter_price(self, ytm_oas_or_ytw):

        if self._clean_price is not None:
            target_value = self._clean_price + self.acc_int(self._as_of_date)
        else:
            target_value = self.dirty_price(self._as_of_date)

        if ytm_oas_or_ytw == "YTM":

            def iter_price(ytm):

                return Bond.__price_lambda(self._as_of_date, "DIRTY", self._maturity_date,
                                           self._coupon_frequency, self._coupon_day_count_conv, self._dated_date,
                                           self._cashflows, self._int_cashflows, self._issue_date, self._notional,
                                           self._amort_schedule, None, self._discount_curve, ytm, self._yield_frequency,
                                           self._yield_day_count_conv, self._first_coupon_date,
                                           self._last_coupon_date, self._reset_frequency, self._ref_rate_tenor,
                                           self._ref_rate_curve, self._ref_rate_frequency,
                                           self._ref_rate_day_count_conv,  self._initial_ref_rate, self._coupon_margin,
                                           None, None, None, self._step_schedule) - target_value

        elif ytm_oas_or_ytw == "OAS":

            def iter_price(oas):

                return Bond.__price_lambda(self._as_of_date, "DIRTY", self._maturity_date,
                                           self._coupon_frequency, self._coupon_day_count_conv, self._dated_date,
                                           self._cashflows, self._int_cashflows, self._issue_date, self._notional,
                                           self._amort_schedule, oas, self._discount_curve, None, self._yield_frequency,
                                           self._yield_day_count_conv, self._first_coupon_date,
                                           self._last_coupon_date, self._reset_frequency, self._ref_rate_tenor,
                                           self._ref_rate_curve, self._ref_rate_frequency,
                                           self._ref_rate_day_count_conv, self._initial_ref_rate, self._coupon_margin,
                                           self._call_schedule, self._put_schedule, None, self._step_schedule) \
                       - target_value

        else:

            def iter_price(ytw):

                return Bond.__price_lambda(self._as_of_date, "DIRTY", self._maturity_date,
                                           self._coupon_frequency, self._coupon_day_count_conv, self._dated_date,
                                           self._cashflows, self._int_cashflows, self._issue_date, self._notional,
                                           self._amort_schedule, None, self._discount_curve, ytw, self._yield_frequency,
                                           self._yield_day_count_conv, self._first_coupon_date,
                                           self._last_coupon_date, self._reset_frequency, self._ref_rate_tenor,
                                           self._ref_rate_curve, self._ref_rate_frequency,
                                           self._ref_rate_day_count_conv, self._initial_ref_rate, self._coupon_margin,
                                           self._call_schedule, self._put_schedule, None, self._step_schedule) \
                       - target_value

        return iter_price

    @staticmethod
    def quick_ir01(val_date, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                   issue_date=None, notional=100, oas=None, discount_curve=None, ytm=None, yield_frequency=None,
                   yield_day_count_conv=None, first_coupon_date=None, last_coupon_date=None, reset_frequency=None,
                   ref_rate_tenor=None, ref_rate_curve=None, ref_rate_frequency=None, ref_rate_day_count_conv=None,
                   initial_ref_rate=None, coupon_margin=None):

        if val_date >= maturity_date:
            return 0.

        dr = 0.0001

        if (oas is not None) and (discount_curve is not None):

            oas -= dr

            if ref_rate_curve is not None:
                ref_rate_curve = ref_rate_curve.parallel_shift(-dr)

            bump_down = Bond.quick_price(val_date, "DIRTY", maturity_date, coupon, coupon_frequency,
                                         coupon_day_count_conv, dated_date, issue_date, notional, oas,
                                         discount_curve, ytm, yield_frequency, yield_day_count_conv,
                                         first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                                         ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate,
                                         coupon_margin)

            oas += 2 * dr

            if ref_rate_curve is not None:
                ref_rate_curve = ref_rate_curve.parallel_shift(2 * dr)

            bump_up = Bond.quick_price(val_date, "DIRTY", maturity_date, coupon, coupon_frequency,
                                       coupon_day_count_conv, dated_date, issue_date, notional, oas,
                                       discount_curve, ytm, yield_frequency, yield_day_count_conv,
                                       first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                                       ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate,
                                       coupon_margin)

        else:

            ytm -= dr

            if ref_rate_curve is not None:
                ref_rate_curve = ref_rate_curve.parallel_shift(-dr)

            bump_down = Bond.quick_price(val_date, "DIRTY", maturity_date, coupon, coupon_frequency,
                                         coupon_day_count_conv, dated_date, issue_date, notional, oas,
                                         discount_curve, ytm, yield_frequency, yield_day_count_conv,
                                         first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                                         ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate,
                                         coupon_margin)

            ytm += 2 * dr

            if ref_rate_curve is not None:
                ref_rate_curve = ref_rate_curve.parallel_shift(2 * dr)

            bump_up = Bond.quick_price(val_date, "DIRTY", maturity_date, coupon, coupon_frequency,
                                       coupon_day_count_conv, dated_date, issue_date, notional, oas,
                                       discount_curve, ytm, yield_frequency, yield_day_count_conv,
                                       first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                                       ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate,
                                       coupon_margin)

        return (bump_down - bump_up) / (2 * dr * 10000)

    def ir01(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        dr = 0.0001

        if (self._oas is not None) and (self._discount_curve is not None):

            self._oas -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)

            self._oas += 2 * dr
            self.__recalculate_cashflows(2 * dr)
            bump_up = self.dirty_price(val_date)

            self._oas -= dr
            self.__recalculate_cashflows(-dr)

        else:

            self._ytm -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)

            self._ytm += 2 * dr
            self.__recalculate_cashflows(2 * dr)
            bump_up = self.dirty_price(val_date)

            self._ytm -= dr
            self.__recalculate_cashflows(-dr)

        return (bump_down - bump_up) / (2 * dr * 10000)

    def __recalculate_cashflows(self, shift):

        if self._ref_rate_curve is not None:

            self._ref_rate_curve = self._ref_rate_curve.parallel_shift(shift)
            self._resets = self._generate_resets()

            self._int_cashflows = self._generate_int_cashflows()
            self._cashflows = self._generate_total_cashflows()

    @staticmethod
    def quick_eff_dur(val_date, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                      issue_date=None, notional=100, oas=None, discount_curve=None, ytm=None, yield_frequency=None,
                      yield_day_count_conv=None, first_coupon_date=None, last_coupon_date=None, reset_frequency=None,
                      ref_rate_tenor=None, ref_rate_curve=None, initial_ref_rate=None, coupon_margin=None):

        if val_date >= maturity_date:
            return 0.

        ir01 = Bond.quick_ir01(val_date, maturity_date, coupon, coupon_frequency, coupon_day_count_conv, dated_date,
                               issue_date, notional, oas, discount_curve, ytm, yield_frequency, yield_day_count_conv,
                               first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor, ref_rate_curve,
                               initial_ref_rate, coupon_margin)

        value = Bond.quick_price(val_date, "DIRTY", maturity_date, coupon, coupon_frequency,
                                 coupon_day_count_conv, dated_date, issue_date, notional, oas,
                                 discount_curve, ytm, yield_frequency, yield_day_count_conv,
                                 first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor, ref_rate_curve,
                                 initial_ref_rate, coupon_margin)

        return (ir01 * 100 / value) * 100

    def eff_dur(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        ir01 = self.ir01(val_date)
        value = self.dirty_price(val_date)

        return (ir01 * 100 / value) * 100

    def ir01_cvx(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        dr = 0.0001
        value = self.dirty_price(val_date)

        if (self._oas is not None) and (self._discount_curve is not None):

            self._oas -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)

            self._oas += 2 * dr
            self.__recalculate_cashflows(2 * dr)
            bump_up = self.dirty_price(val_date)

            self._oas -= dr
            self.__recalculate_cashflows(-dr)

        else:

            self._ytm -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)

            self._ytm += 2 * dr
            self.__recalculate_cashflows(2 * dr)
            bump_up = self.dirty_price(val_date)

            self._ytm -= dr
            self.__recalculate_cashflows(-dr)

        return (bump_up - (2 * value) + bump_down) / (dr * 10000)**2

    def eff_cvx(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        ir01_cvx = self.ir01_cvx(val_date)
        value = self.dirty_price(val_date)

        return (ir01_cvx * 10000 / value) * 100

    def mod_dur(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        dr = 0.01

        value = self.dirty_price(val_date)

        if (self._oas is not None) and (self._discount_curve is not None):

            self._oas -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)
            self._oas += dr
            self.__recalculate_cashflows(dr)

        else:

            self._ytm -= dr
            self.__recalculate_cashflows(-dr)
            bump_down = self.dirty_price(val_date)
            self._ytm += dr
            self.__recalculate_cashflows(dr)

        return ((bump_down - value) / value) * 100

    def macaulay_dur(self, val_date):

        if val_date >= self._maturity_date:
            return 0.

        if self._ytm is not None:

            dur = 0.
            value = self.dirty_price(val_date)

            df = self.__ytm_df(val_date, self._ytm, self._yield_day_count_conv, self._yield_frequency,
                               self._maturity_date)

            for cf_date, cf_val in self._cashflows.items():

                if cf_date > val_date:

                    t = Date.year_frac(val_date, cf_date, self._coupon_day_count_conv, self._maturity_date)
                    dur += (t * cf_val * df(cf_date)) / value

            return dur

        else:
            raise ValueError("YTM must be known to calculate Macaulay duration")


class FixedRateBond(Bond):
    """Class for fixed rate bonds"""

    def __init__(self, name, maturity_date, coupon, frequency, day_count_conv, dated_date, issue_date=None,
                 notional=100., clean_price=None, as_of_date=None, amort_schedule=None, oas=None, discount_curve=None,
                 ytm=None, yield_frequency=None, yield_day_count_conv=None, first_coupon_date=None,
                 last_coupon_date=None, call_schedule=None, put_schedule=None, ytw=None, step_schedule=None):

        super().__init__(name, maturity_date, coupon, frequency, day_count_conv, dated_date, issue_date, notional,
                         clean_price, as_of_date, amort_schedule, oas, discount_curve, ytm, yield_frequency,
                         yield_day_count_conv, first_coupon_date, last_coupon_date, None, None, None, None, None, None,
                         None, call_schedule, put_schedule, ytw, step_schedule, None)

    @staticmethod
    def _calc_coupon_lambda(date, notional=None, coupon=None, step_schedule=None):

        if coupon is None:
            return 0.

        if step_schedule is not None:

            dates = list(step_schedule.keys())
            rates = list(step_schedule.values())

            for i in range(0, len(dates)):

                if date <= dates[i]:
                    return notional * (rates[i] if (i == 0) else rates[i - 1])

        return notional * coupon

    def _generate_int_cashflows(self):

        int_cashflows = Bond._generate_int_cashflows_lambda(self._int_cashflow_schedule, self._prin_schedule,
                                                            self._coupon, self._resets, self._coupon_margin,
                                                            self._eff_first_coupon_date, self._eff_last_coupon_date,
                                                            self._coupon_day_count_conv,
                                                            FixedRateBond._calc_coupon_lambda, self._step_schedule,
                                                            None)

        self._int_cashflow_values = int_cashflows.values()
        return int_cashflows


class FloatingRateBond(Bond):
    """Class for floating rate bonds"""

    def __init__(self, name, maturity_date, frequency, day_count_conv, reset_frequency, ref_rate_tenor, ref_rate_curve,
                 ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate, coupon_margin, dated_date,
                 issue_date=None, notional=100., clean_price=None, as_of_date=None, amort_schedule=None, oas=None,
                 discount_curve=None, ytm=None, yield_frequency=None, yield_day_count_conv=None, first_coupon_date=None,
                 last_coupon_date=None, call_schedule=None, put_schedule=None, ytw=None):

        super().__init__(name, maturity_date, None, frequency, day_count_conv, dated_date, issue_date, notional,
                         clean_price, as_of_date, amort_schedule, oas, discount_curve, ytm, yield_frequency,
                         yield_day_count_conv, first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                         ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate, coupon_margin,
                         call_schedule, put_schedule, ytw, None, None)

    @staticmethod
    def _calc_coupon_lambda(date, notional=None, resets=None, coupon_margin=None):

        dates = list(resets.keys())
        rates = list(resets.values())

        for i in range(0, len(dates)):

            if date <= dates[i]:
                return notional * ((coupon_margin + rates[i]) if (i == 0) else (coupon_margin + rates[i - 1]))

        return 0.

    def _generate_int_cashflows(self):

        int_cashflows = Bond._generate_int_cashflows_lambda(self._int_cashflow_schedule, self._prin_schedule,
                                                            self._coupon, self._resets, self._coupon_margin,
                                                            self._eff_first_coupon_date, self._eff_last_coupon_date,
                                                            self._coupon_day_count_conv,
                                                            FloatingRateBond._calc_coupon_lambda, self._step_schedule,
                                                            None)

        self._int_cashflow_values = int_cashflows.values()

        return int_cashflows


class FixedToFloatingRateBond(Bond):
    """Class for fixed-to-floating rate bonds"""

    def __init__(self, name, maturity_date, coupon, frequency, day_count_conv, reset_frequency, ref_rate_tenor,
                 ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, coupon_margin, first_floating_date,
                 dated_date, issue_date=None, notional=100., clean_price=None, as_of_date=None, amort_schedule=None,
                 initial_ref_rate=None, oas=None, discount_curve=None, ytm=None, yield_frequency=None,
                 yield_day_count_conv=None, first_coupon_date=None, last_coupon_date=None, call_schedule=None,
                 put_schedule=None, ytw=None, step_schedule=None):

        super().__init__(name, maturity_date, coupon, frequency, day_count_conv, dated_date, issue_date, notional,
                         clean_price, as_of_date, amort_schedule, oas, discount_curve, ytm, yield_frequency,
                         yield_day_count_conv, first_coupon_date, last_coupon_date, reset_frequency, ref_rate_tenor,
                         ref_rate_curve, ref_rate_frequency, ref_rate_day_count_conv, initial_ref_rate, coupon_margin,
                         call_schedule, put_schedule, ytw, step_schedule, first_floating_date)

    @staticmethod
    def _calc_coupon_lambda(date, notional=None, coupon=None, step_schedule=None, resets=None, coupon_margin=None,
                            first_floating_date=None):

        dates = list(resets.keys())
        rates = list(resets.values())

        if date < first_floating_date:

            if step_schedule is not None:

                dates = list(step_schedule.keys())
                rates = list(step_schedule.values())

                for i in range(0, len(dates)):

                    if date <= dates[i]:
                        return notional * (rates[i] if (i == 0) else rates[i - 1])

            return notional * coupon
        
        for i in range(0, len(dates)):

            if date <= dates[i]:
                return notional * ((coupon_margin + rates[i]) if (i == 0) else (coupon_margin + rates[i - 1]))

        return 0.

    def _generate_int_cashflows(self):

        int_cashflows = Bond._generate_int_cashflows_lambda(self._int_cashflow_schedule, self._prin_schedule,
                                                            self._coupon, self._resets, self._coupon_margin,
                                                            self._eff_first_coupon_date, self._eff_last_coupon_date,
                                                            self._coupon_day_count_conv,
                                                            FixedToFloatingRateBond._calc_coupon_lambda,
                                                            self._step_schedule, self._first_floating_date)

        self._int_cashflow_values = int_cashflows.values()

        return int_cashflows

'''
maturity = datetime.date(2020, 8, 24)
dated_date = datetime.date(2010, 8, 24)
issue_date = datetime.date(2010, 8, 24)
first_coupon_date = datetime.date(2016, 1, 31)
last_coupon_date = datetime.date(2025, 2, 28)

test_bond = FixedRateBond("test_bond", maturity, .04875, "S", "30/360", dated_date, None, 100., None, None, None, None,
                          None, "S", "30/360", None, None)
cashflows = test_bond.cashflows

print(Date.year_frac(datetime.date(2016, 6, 30), datetime.date(2020, 8, 24), "30/360"))

for key, value in cashflows.items():
    print(str(key) + ": " + str(value))
'''
