#!/bin/bash
source bin/functions.sh

#? some statements here that will probably be used
PATHING_SINGLE="FALSE"
PATHING_DIFFERENT="FALSE"
db_path_single_confirmed="no"


# ! start the actual script here

database_installer_intro

echo -e "Downloading and installing the various databases will take several hours depending on your local setup"
while read -r -p "Do you wish to continue with this installation process now? [y/N] " db_continue_response
do 
    db_continue_response=${db_continue_response,,}
    if [[ "${db_continue_response}" =~ ^(no|n)$ ]]; then
        minispacer
        echo -e "Exiting the Database installation process on user request..."
        sleep 1
        exit 0
    elif [[ "${db_continue_response}" =~ ^(yes|y)$ ]]; then
        minispacer
        echo -e "Continuing with the installation process..."
        sleep 1
        break
    else
        echo -e "Please answer with 'yes' or 'no'"
    fi
done

database_installer

# ! Introduction to paths for the database installer
echo -e "Do you wish to use \e[1mone single path for all databases\e[0m or do you wish to have \e[1mindividual paths for each database?\e[0m"
thinline
echo -e "Example for one single path for all databases:"
echo -e "/same/path/to/database_A/"
echo -e "/same/path/to/database_B/"
echo -e "/same/path/to/database_C/"
thinline
echo -e "Examples for individual paths for each database:"
echo -e "/first/path/to/database_A/"
echo -e "/different/path/to/other/database_B/"
echo -e "/third/path/to/database_C/"
minispacer



while read -r -p "Which path-structure do you prefer for your system? Single or Individual? [S/I] " db_pathstructure_response
do
    db_pathstructure_response=${db_pathstructure_response,,}
    if [[ "${db_pathstructure_response}" =~ ^(single|s)$ ]]; then
        PATHING_SINGLE="TRUE"
        echo "A single path-structure has been selected."
        echo "Continuing..."
        sleep 2
        break
    elif [[ "${db_pathstructure_response}" =~ ^(individual|i)$ ]]; then
        PATHING_DIFFERENT="TRUE"
        echo "An individual path-structure for each database has been selected."
        echo "Continuing..."
        sleep 2
        break
    else
        echo -e "Please answer with 'single' or 'individual'"
    fi
done


database_installer
if [ "${PATHING_SINGLE}" == "TRUE" ]; then
    # ! ASK THE PROMPT FOR THE BASE PATH HERE
    echo -e "The databases have to be installed in a location from the root directory..."
    echo -e "In other words, the path you specify must start with a '/'"
    echo -e "It is also important to note that the specified path should NOT END with a trailing '/'"
    thinline
    echo -e "Example of a path which is correct:"
    echo -e "\e[40;92m /this/path/is/correct \e[0m"
    thinline
    echo -e "Example of a path which is incorrect:"
    echo -e "\e[40;91m this/path/is/incorrect/ \e[0m"
    minispacer

    while read -r -p "Please type out the path under which you want to install the databases: " db_path_single_response
    do 
        db_path_single_response=${db_path_single_response}
        minispacer
        echo -e "With the path you specified, the database locations will look like this:"
        echo -e "\e[1m${db_path_single_response}/NT_database\e[0m"
        echo -e "\e[1m${db_path_single_response}/NR_database\e[0m"
        echo -e "\e[1m${db_path_single_response}/taxdb\e[0m"
        echo -e "\e[1m${db_path_single_response}/new_taxdump\e[0m"
        echo -e "\e[1m${db_path_single_response}/krona_taxonomy\e[0m"
        echo -e "\e[1m${db_path_single_response}/Virus-Host_interaction_db\e[0m"
        minispacer
        
        while read -r -p "Does this look correct? [yes/no] " db_path_single_confirmed
        do
            db_path_single_confirmed=${db_path_single_confirmed,,}
            if [[ "${db_path_single_confirmed}" =~ ^(yes)$ ]]; then
                echo "Paths for the databases have been confirmed"
                echo "Continuing..."
                sleep 5
                break
            elif [[ "${db_path_single_confirmed}" =~ ^(no)$ ]]; then
                minispacer
                break
            else
                echo "Please answer with 'yes' or 'no'"
            fi
        done
        
        if [ "${db_path_single_confirmed}" == "yes" ]; then
        
        database_installer

        # ! Set paths
        DB_PATH_NT="${db_path_single_response}/NT_database"
        DB_PATH_NR="${db_path_single_response}/NR_database"
        DB_PATH_TAX="${db_path_single_response}/taxdb"
        DB_PATH_NTAX="${db_path_single_response}/new_taxdump"
        DB_PATH_KRONA="${db_path_single_response}/krona_taxonomy"
        DB_PATH_VHOST="${db_path_single_response}/Virus-Host_interaction_db"
        
        # ! Make folders
        mkdir -p ${DB_PATH_NT}
        mkdir -p ${DB_PATH_NR}
        mkdir -p ${DB_PATH_TAX}
        mkdir -p ${DB_PATH_NTAX}
        mkdir -p ${DB_PATH_KRONA}
        mkdir -p ${DB_PATH_VHOST}
        
        # ! Start downloading databases
            # ! download taxdb
            cd "${DB_PATH_TAX}"
            printf "\nDownloading taxonomy database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress taxdb
            
            # ! download krona LCA db
            cd "${DB_PATH_KRONA}"
            printf "\nDownloading Krona LCA database... \n"
            bash "${CONDA_PREFIX}"/opt/krona/updateTaxonomy.sh ./
            bash "${CONDA_PREFIX}"/opt/krona/updateAccessions.sh ./

            # ! download VHOST-interaction db
            cd "${DB_PATH_VHOST}"
            printf "\nDownloading Virus-Host interaction database... \n"
            curl -o virushostdb.tsv -L ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv            
            
            # ! download new_taxdump db
            cd "${DB_PATH_NTAX}"
            printf "\nDownloading NCBI New_taxdump database... \n"
            curl -o new_taxdump.tar.gz -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
            curl -o new_taxdump.tar.gz.md5 -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
            tar -xzf new_taxdump.tar.gz
            for file in *.dmp;
                do awk '{gsub("\t",""); if(substr($0,length($0),length($0))=="|") print substr($0,0,length($0)-1); else print $0}' < ${file} > ${file}.delim;
                done
                
            # ! download BLAST NT db
            cd "${DB_PATH_NT}"
            printf "\nDownloading NCBI BLAST NT database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nt
            
            # ! download BLAST NR db
            cd "${DB_PATH_NR}"
            printf "\nDownloding NCBI BLAST NR database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nr

            sleep 2
            printf "\nWe donezo \n"
            printf "\n :) \n"

        break
        fi

    done
fi

if [ "${PATHING_DIFFERENT}" == "TRUE" ]; then
    # ! ASK THE DIFFERENT PROMPTS FOR EACH DATABASE HERE
    echo "this part isn't done yet"
    echo ":)"
fi



#//     echo -e "${db_path_single_response}"

#//     echo -e "${PATHING_SINGLE}"
#//     echo -e "${PATHING_DIFFERENT}"












#//if [[ "${db_continue_response}" =~ ^(no|n)$ ]]; then
#//    minispacer
#//    echo -e "Exiting the Database installation process on user request..."
#//    sleep 1
#//    exit 0
#//elif [[ "${db_continue_response}" =~ ^(yes|y)$ ]]; then
#//    minispacer
#//    echo -e "Continuing with the installation process..."
#//    sleep 1
#//else
#//    read -r -p "Please answer with 'yes' or 'no' [y/N] " db_continue_response
#//
#//fi