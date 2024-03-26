wget -c wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 -O resources/strelka-2.9.10.centos6_x86_64.tar.bz2
cd resources
tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2

#Test installation
bash strelka-2.9.10.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.10.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash

cd -
