.PHONY: make_conda update_conda remove_conda update_data update_data_dry update_data_dry

.ONESHELL:

PROJECT_NAME = string_sims

make_conda:
	conda env create -f environment.yml
	conda activate $(PROJECT_NAME)
	ipython kernel install --user --name=$(PROJECT_NAME)

update_conda:
	conda env update --file environment.yml
	ipython kernel install --user --name=$(PROJECT_NAME)

remove_conda:
	conda remove --name=$(PROJECT_NAME) --all

update_data:
	 while read a b;do rsync -rauLih --progress --delete --exclude-from=data/raw/exclude_list.txt   $$a $$b;done < data/raw/dir_list.txt
	 python  src/data/get_external_data.py

update_data_dry:
	 while read a b;do rsync -raunLih --delete --exclude-from=data/raw/exclude_list.txt   $$a $$b;done < data/raw/dir_list.txt

update_data_dry_verbose:
	 while read a b;do rsync -raunLivvvh --delete --exclude-from=data/raw/exclude_list.txt   $$a $$b;done < data/raw/dir_list.txt

clean:
	find . -iname \#* -exec rm {} \;
	find . -iname "slurm*err" -exec rm {} \;
	find . -iname "slurm*out" -exec rm {} \;

format:
	black -l 79 .
	isort .

help:
	@echo "Possible options:"
	@echo "make_conda"
	@echo "update_conda"
	@echo "remove_conda"
	@echo "update_data"
	@echo "update_data_dry"
