token=$(<gdc-user-token.2020-11-04T23_18_05.325Z.txt)

curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/e1e10446-ae53-45b5-8160-d1dd40637cff
?gencode=TERT&gencode=TRIO' --output get_examp_slice.bam
