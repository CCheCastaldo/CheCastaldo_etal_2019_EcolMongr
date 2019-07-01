#!/bin/bash

# nohup bash /home/JobRemote.bash > foo.out 2> foo.err < /dev/null &

run=GlobalModel

# initialize git
git config --global user.name "CCheCastaldo"
git config --global user.email "cccweb@icloud.com"
git config --global push.default simple
 echo 'Host github.com
	StrictHostKeyChecking no' > ~/.ssh/config


# clone repo
cd /home
git clone git@github.com:CCheCastaldo/SalixStressVigor.git

 
# run R job
cd /home/SalixStressVigor/ModelBuild/$run
model_id=$(wget -q -O - http://instance-data/latest/meta-data/instance-id)
Rscript --no-save --no-restore --verbose Rsetup.R
Rscript --no-save --no-restore --verbose Implementation.R

# put run results in subfolder with instance_id as name
mkdir $model_id
cp *.Rmd $model_id
cp Implementation.R $model_id
mv *.rda $model_id

# commit changes to repo
git add --all
git commit -m "ec2 $run $model_id run complete"
git pull
git push origin

# self destruct
#sudo shutdown -h now

