#!/bin/bash
# shellcheck disable=SC1091

#* load external functions
source bin/functions.sh


if [ "${1}" == "master" ]; then
	echo -e "Updating Jovian will delete \e[1m all local changes\e[0m that you've made.\nIf you have any important changes then please make a backup of the files you changed before you continue"
	while read -r -p "Are you sure you wish to update Jovian now? [y/N] " updateresponse
	do
		updateresponse=${updateresponse,,}
		if [[ "${updateresponse}" =~ ^(yes|y)$ ]]; then
			break
		elif [[ "${updateresponse}" =~ ^(no|n)$ ]]; then
			echo -e "Aborting update process on user request"
			exit 0
		else
			echo -e "Please answer with 'yes' or 'no'"
		fi
	done
		
	git reset --hard
	git checkout master

	echo -e "DONE"
else

	echo -e "Changing the Jovian version will delete \e[1m all local changes\e[0m that you've made. If you have any important changes then please make a backup of the files you changed before you continue"
	while read -r -p "Are you sure you wish to change the Jovian version? [y/N] " changeresponse
	do
		changeresponse=${changeresponse,,}
		if [[ "${changeresponse}" =~ ^(yes|y)$ ]]; then
			break
		elif [[ "${changeresponse}" =~ ^(no|n)$ ]]; then
			echo -e "Aborting the version change process on user request"
			exit 0
		else
			echo -e "Please answer with 'yes' or 'no'"
		fi
	done
	git reset --hard
	git checkout master	
	git fetch origin --tags "${1}"
	git reset --hard FETCH_HEAD

	echo -e "DONE"
fi