#!/bin/bash
# shellcheck disable=SC1091
source bin/functions.sh

#? some statements here that will probably be used
PATHING_SINGLE="FALSE"
PATHING_DIFFERENT="FALSE"
db_path_single_confirmed="no"
db_paths_indiv_confirmed="no"


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

    while read -r -e -p "Please type out the path under which you want to install the databases: " db_path_single_response
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
            cd "${DB_PATH_TAX}" || exit
            printf "\nDownloading taxonomy database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress taxdb
            
            # ! download krona LCA db
            cd "${DB_PATH_KRONA}" || exit
            printf "\nDownloading Krona LCA database... \n"
            bash "${CONDA_PREFIX}"/opt/krona/updateTaxonomy.sh ./
            bash "${CONDA_PREFIX}"/opt/krona/updateAccessions.sh ./

            # ! download VHOST-interaction db
            cd "${DB_PATH_VHOST}" || exit
            printf "\nDownloading Virus-Host interaction database... \n"
            curl -o virushostdb.tsv -L ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv            
            
            # ! download new_taxdump db
            cd "${DB_PATH_NTAX}" || exit
            printf "\nDownloading NCBI New_taxdump database... \n"
            curl -o new_taxdump.tar.gz -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
            curl -o new_taxdump.tar.gz.md5 -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
            tar -xzf new_taxdump.tar.gz
            for file in *.dmp;
                do awk '{gsub("\t",""); if(substr($0,length($0),length($0))=="|") print substr($0,0,length($0)-1); else print $0}' < ${file} > ${file}.delim;
                done
                
            # ! download BLAST NT db
            cd "${DB_PATH_NT}" || exit
            printf "\nDownloading NCBI BLAST NT database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nt
            
            # ! download BLAST NR db
            cd "${DB_PATH_NR}" || exit
            printf "\nDownloding NCBI BLAST NR database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nr

            sleep 2
            printf "\nDone with downloading and installing the databases required for Jovian \n"
            printf "The databases have been placed in the %s folder\n" "${db_path_single_response}"

            cat << EOF > ~/database-updater.sh
#!/bin/bash
source activate Jovian_helper

### Updating BLAST taxdb
cd "${DB_PATH_TAX}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress taxdb

### UPDATING KRONA
cd "${DB_PATH_KRONA}" || exit 1
bash "${CONDA_PREFIX}"/opt/krona/updateTaxonomy.sh ./
bash "${CONDA_PREFIX}"/opt/krona/updateAccessions.sh ./

### UPDATING VHOST INTERACTION DB
cd "${DB_PATH_VHOST}" || exit 1
curl -o virushostdb.tsv -L ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv

### UPDATING NEWTAXDUMP DB
cd "${DB_PATH_NTAX}" || exit 1
curl -o new_taxdump.tar.gz -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
curl -o new_taxdump.tar.gz.md5 -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
tar -xzf new_taxdump.tar.gz
    for file in *.dmp;
        do awk '{gsub("\t",""); if(substr($0,length($0),length($0))=="|") print substr($0,0,length($0)-1); else print $0}' < ${file} > ${file}.delim;
        done

### UPDATING BLAST NT
cd "${DB_PATH_NT}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nt

### UPDATING BLAST NR
cd "${DB_PATH_NR}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nr

echo -e "Done with updating databases"
exit

EOF
        
        break
        fi

    done
fi

if [ "${PATHING_DIFFERENT}" == "TRUE" ]; then
    # ! ASK THE DIFFERENT PROMPTS FOR EACH DATABASE HERE
    printf "\nThe databases have to be installed in a location from the root directory...\n"
    printf "In other words, the paths you specify must start with a '/'\n"
    printf "It is also important to note that the specified paths should NOT END with a trailing '/'\n"
    thinline
    printf "Example of paths which are correct:\n"
    printf "\e[40;92m /this/path/is/correct \e[0m\n"
    printf "\e[40;92m /this/path/is/also/correct \e[0m\n"
    thinline
    printf "Example of paths which are incorrect:\n"
    printf "\e[40;91m this/path/is/incorrect/ \e[0m\n"
    printf "\e[40;91m this/path/is/also/incorrect \e[0m\n"
    minispacer

    while read -r -e -p "Please type out the path under which you want to install the BLAST NT database: " db_path_indiv_NT_response; read -r -e -p "Please type out the path under which you want to install the BLAST NR database: " db_path_indiv_NR_response; read -r -e -p "Please type out the path under which you want to install the TaxDB: " db_path_indiv_taxdb_response; read -r -e -p "Please type out the path under which you want to install the NewTaxdump database: " db_path_indiv_ntax_response; read -r -e -p "Please type out the path under which you want to install the Krona LCA database: " db_path_indiv_krona_response; read -r -e -p "please type out the path under which you want to install the VirusHost Interaction database: " db_path_indiv_vhost_response
    do
        db_path_indiv_NT_response=${db_path_indiv_NT_response}
        db_path_indiv_NR_response=${db_path_indiv_NR_response}
        db_path_indiv_taxdb_response=${db_path_indiv_taxdb_response}
        db_path_indiv_ntax_response=${db_path_indiv_ntax_response}
        db_path_indiv_krona_response=${db_path_indiv_krona_response}
        db_path_indiv_vhost_response=${db_path_indiv_vhost_response}

        minispacer
        printf "\nWith the paths you specified, the database locations will look like this:\n"
        printf "\e[1m%s/NT_database\e[0m\n" "${db_path_indiv_NT_response}"
        printf "\e[1m%s/NR_database\e[0m\n" "${db_path_indiv_NR_response}"
        printf "\e[1m%s/taxdb\e[0m\n" "${db_path_indiv_taxdb_response}"
        printf "\e[1m%s/new_taxdump\e[0m\n" "${db_path_indiv_ntax_response}"
        printf "\e[1m%s/krona_taxonomy\e[0m\n" "${db_path_indiv_krona_response}"
        printf "\e[1m%s/Virus-Host_interaction_db\e[0m\n" "${db_path_indiv_vhost_response}"
        minispacer

        while read -r -p "Does this look correct? [yes/no] " db_paths_indiv_confirmed
        do
            db_paths_indiv_confirmed=${db_paths_indiv_confirmed,,}
            if [[ "${db_paths_indiv_confirmed}" =~ ^(yes)$ ]]; then
                echo "Paths for the databases have been confirmed"
                echo "Continuing..."
                sleep 5
                break
            elif [[ "${db_paths_indiv_confirmed}" =~ ^(no)$ ]]; then
                minispacer
                break
            else
                echo "Please answer with 'yes' or 'no'"
            fi
        done
        
        if [ "${db_paths_indiv_confirmed}" == "yes" ]; then

        database_installer
            
        # ! Set paths
        DB_PATH_NT="${db_path_indiv_NT_response}/NT_database"
        DB_PATH_NR="${db_path_indiv_NR_response}/NR_database"
        DB_PATH_TAX="${db_path_indiv_taxdb_response}/taxdb"
        DB_PATH_NTAX="${db_path_indiv_ntax_response}/new_taxdump"
        DB_PATH_KRONA="${db_path_indiv_krona_response}/krona_taxonomy"
        DB_PATH_VHOST="${db_path_indiv_vhost_response}/Virus-Host_interaction_db"


        # ! Make folders
        mkdir -p ${DB_PATH_NT}
        mkdir -p ${DB_PATH_NR}
        mkdir -p ${DB_PATH_TAX}
        mkdir -p ${DB_PATH_NTAX}
        mkdir -p ${DB_PATH_KRONA}
        mkdir -p ${DB_PATH_VHOST}

        # ! Start downloading databases
            # ! Download taxdb
            cd "${DB_PATH_TAX}" || exit
            printf "\nDownloading taxonomy database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress taxdb
                
            # ! download krona LCA db
            cd "${DB_PATH_KRONA}" || exit
            printf "\nDownloading Krona LCA database... \n"
            bash "${CONDA_PREFIX}"/opt/krona/updateTaxonomy.sh ./
            bash "${CONDA_PREFIX}"/opt/krona/updateAccessions.sh ./

            # ! download VHOST-interaction db
            cd "${DB_PATH_VHOST}" || exit
            printf "\nDownloading Virus-Host interaction database... \n"
            curl -o virushostdb.tsv -L ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv            
                
            # ! download new_taxdump db
            cd "${DB_PATH_NTAX}" || exit
            printf "\nDownloading NCBI New_taxdump database... \n"
            curl -o new_taxdump.tar.gz -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
            curl -o new_taxdump.tar.gz.md5 -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
            tar -xzf new_taxdump.tar.gz
            for file in *.dmp;
                do awk '{gsub("\t",""); if(substr($0,length($0),length($0))=="|") print substr($0,0,length($0)-1); else print $0}' < ${file} > ${file}.delim;
                done
                    
            # ! download BLAST NT db
            cd "${DB_PATH_NT}" || exit
            printf "\nDownloading NCBI BLAST NT database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nt
                
            # ! download BLAST NR db
            cd "${DB_PATH_NR}" || exit
            printf "\nDownloding NCBI BLAST NR database... \n"
            perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nr

            sleep 2
            printf "\nDone with downloading and installing the databases required for Jovian \n"
            printf "The databases have been placed in the following folders:\n"
            printf "\e[1m%s/NT_database\e[0m\n" "${db_path_indiv_NT_response}"
            printf "\e[1m%s/NR_database\e[0m\n" "${db_path_indiv_NR_response}"
            printf "\e[1m%s/taxdb\e[0m\n" "${db_path_indiv_taxdb_response}"
            printf "\e[1m%s/new_taxdump\e[0m\n" "${db_path_indiv_ntax_response}"
            printf "\e[1m%s/krona_taxonomy\e[0m\n" "${db_path_indiv_krona_response}"
            printf "\e[1m%s/Virus-Host_interaction_db\e[0m\n" "${db_path_indiv_vhost_response}"

            cat << EOF > ~/database-updater.sh
#!/bin/bash
source activate Jovian_helper

### Updating BLAST taxdb
cd "${DB_PATH_TAX}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress taxdb

### UPDATING KRONA
cd "${DB_PATH_KRONA}" || exit 1
bash "${CONDA_PREFIX}"/opt/krona/updateTaxonomy.sh ./
bash "${CONDA_PREFIX}"/opt/krona/updateAccessions.sh ./

### UPDATING VHOST INTERACTION DB
cd "${DB_PATH_VHOST}" || exit 1
curl -o virushostdb.tsv -L ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv

### UPDATING NEWTAXDUMP DB
cd "${DB_PATH_NTAX}" || exit 1
curl -o new_taxdump.tar.gz -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
curl -o new_taxdump.tar.gz.md5 -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
tar -xzf new_taxdump.tar.gz
    for file in *.dmp;
        do awk '{gsub("\t",""); if(substr($0,length($0),length($0))=="|") print substr($0,0,length($0)-1); else print $0}' < ${file} > ${file}.delim;
        done

### UPDATING BLAST NT
cd "${DB_PATH_NT}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nt

### UPDATING BLAST NR
cd "${DB_PATH_NR}" || exit 1
perl "${CONDA_PREFIX}"/bin/update_blastdb.pl --decompress nr

echo -e "Done with updating databases"
exit

EOF

        break
        fi

    done
fi

