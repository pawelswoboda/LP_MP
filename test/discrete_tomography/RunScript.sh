# This runs DT, so that you can see the message issue

#
# in see discrete_tomography_message_counting_pairwise.hxx line 111 we output
# every reparametrization. The assumption is that every second line should only contain
# 0 or inf.
#

./discrete_tomography_debug --tightenConstraintsMax 1 -i GM_low_sp.txt --protocolateConsole | tee log.txt
