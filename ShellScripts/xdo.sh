#!/bin/bash
#takes PhoneNumber as $1 and activates skype window with that name
#then activates dialpad and punches the sequence of dialpad numbers

#----------------------------------------------------------------------------------------
# old school "skype-dialer"
#----------------------------------------------------------------------------------------

#DETECT SOUND WITH PYTHON, WAIT 3s, THEN RUN THIS

PhoneNumber=$1
IdNumber=$2

#Coordinates of Skype Call window:
activate_dialpad=$'350 600'
hangup=$'464 600'
one=$'300 420'
away=$'350 800'

#declare array of dialpad numbers: 0 1 2 3 4 5 6 7 8 9 
declare -a dialpad=("350 550" "300 420" "350 420" "405 420" "300 470" "350 470" "405 470" "300 500" "350 500" "405 500")

#move mouse to "away" location
xdotool mousemove --sync $away 
sleep 0.5

#activate "Call with ..." window
echo "Activating Call window"
echo "Call with +$PhoneNumber"
xdotool search --name "Call with +$PhoneNumber" windowactivate --sync windowsize --sync 557 574 windowmove --sync 0 0
sleep 0.5s


#activate dialpad
#echo "Activating dialpad"
xdotool mousemove --sync $activate_dialpad click 1
sleep 0.5s

#press one for English
xdotool mousemove --sync $one click 1

#wait after "Enter your Id Number"
sleep 3s

#enter IdNumber on the dialpad

#itterate over the digits in IdNumber
for ((i=0; i<${#IdNumber}; i++)) do
   dig=${IdNumber:$i:1} #digit
   #access coordinates of the digit
   xy=${dialpad[$dig]} 
   #echo "dig= $dig xy= $xy"
   #press digit on the dialpad
   xdotool mousemove --sync $xy click 1
   sleep 0.5s 
done

# wait after "Does your last name begin with ..."
sleep 15s
# press one
xdotool mousemove --sync $one click 1

#move mouse to "away" location
xdotool mousemove --sync $away 


#RECORD WITH PYTHON...

#SCRATCH----------------------------------------------------------------------------

#access one element of an array (starts from 0'th index)
#echo "${dialpad[0]}"


#for i in "${dialpad[@]}"; do
#   echo "$i"
#done

