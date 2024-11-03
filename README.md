# EEL3111C
This repo contains the libraries I've created for the TI-NSpire for EEL3111C exams.

**For Exam 2** _updated 11/03 @ 14:23_<br>
**To run, open file, go to page 3, and press "CTRL+R"**<br><br>
Contents
- rlc_response()
  > Inputs: a circuit's resistance, capacitance, and inducatance and if it is parallel/series and if looking for natural/step response<br>
  > Outputs: whether overdamped, underdamped, critically damped, roots of char eq, neper frequency, resonant frequency, all associated equations

- Notes for filter-related problems
- Miscellaneous notes

**Important:**<br>
The TI-NSpire runs an older version of python, and therefore does not support the following python elements
> Match/Case statements (use if/elif/else instead)<br>
> fstrings (alternative example: myStr = "value1 = {}, value2 = {}".format(v1, v2) )<br>
> `__name__ == __main__` (just put code to be automatically run at the very end of the .py file outside any functions)
