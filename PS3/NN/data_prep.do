import delimited "C:\Users\vaguiar\Downloads\TopicsDecisionMaking\NN\data\ABKK_full_experiment.csv", encoding(UTF-8) clear
drop filter1 notus repeatedexcludingotherfilters repeatedandsamedemographics repetitions2 ip

//gen attributes
gen prob50=0
gen prob48=0
gen prob30=0
gen prob14=0
gen prob12=0
gen prob10=0
gen prob0=0
// replace 
replace prob50=1/2 if choice==1
replace prob50=1/4 if choice==3 | choice==4
replace prob48=1/5 if choice==4 | choice==5
replace prob30=1/2 if choice==2
replace prob30=1/4 if choice==3 | choice==5
replace prob14=3/20 if choice==4 | choice==5
replace prob12=1 if choice==0
replace prob10=1/2 if choice==2
replace prob10=1/4 if choice==3 | choice==5
replace prob0=1/2 if choice==1
replace prob0=1/4 if choice==3
replace prob0=2/5 if choice==4
replace prob0=3/20 if choice==5

//Expectations
gen expect=0
replace expect=25 if choice==1
replace expect=20 if choice==2
replace expect=22.5 if choice==3
replace expect=24.125 if choice==4
replace expect=21.625 if choice==5
replace expect=12 if choice==0
//Complexity
gen complex=0
replace complex=1 if choice==1
replace complex=1 if choice==2
replace complex=3 if choice==3
replace complex=3 if choice==4
replace complex=4 if choice==5
replace complex=0 if choice==0
//Position


//export
export delimited using "C:\Users\vaguiar\Downloads\TopicsDecisionMaking\NN\data\ABKK_attributes.csv", replace






