zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-EnR_BETA,ExR-EnR_SE,ExR-EnR_P | sed -e 's/ExR-EnR_//g' | gzip -c > XiongZ_31763980.ExR-EnR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-EnL_BETA,ExR-EnL_SE,ExR-EnL_P | sed -e 's/ExR-EnL_//g' | gzip -c > XiongZ_31763980.ExR-EnL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-ExL_BETA,ExR-ExL_SE,ExR-ExL_P | sed -e 's/ExR-ExL_//g' | gzip -c > XiongZ_31763980.ExR-ExL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-EnL_BETA,EnR-EnL_SE,EnR-EnL_P | sed -e 's/EnR-EnL_//g' | gzip -c > XiongZ_31763980.EnR-EnL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-ExL_BETA,EnR-ExL_SE,EnR-ExL_P | sed -e 's/EnR-ExL_//g' | gzip -c > XiongZ_31763980.EnR-ExL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-ExL_BETA,EnL-ExL_SE,EnL-ExL_P | sed -e 's/EnL-ExL_//g' | gzip -c > XiongZ_31763980.EnL-ExL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-Prn_BETA,N-Prn_SE,N-Prn_P | sed -e 's/N-Prn_//g' | gzip -c > XiongZ_31763980.N-Prn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-AlL_BETA,N-AlL_SE,N-AlL_P | sed -e 's/N-AlL_//g' | gzip -c > XiongZ_31763980.N-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-Sn_BETA,N-Sn_SE,N-Sn_P | sed -e 's/N-Sn_//g' | gzip -c > XiongZ_31763980.N-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-AlR_BETA,N-AlR_SE,N-AlR_P | sed -e 's/N-AlR_//g' | gzip -c > XiongZ_31763980.N-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-AlL_BETA,Prn-AlL_SE,Prn-AlL_P | sed -e 's/Prn-AlL_//g' | gzip -c > XiongZ_31763980.Prn-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-Sn_BETA,Prn-Sn_SE,Prn-Sn_P | sed -e 's/Prn-Sn_//g' | gzip -c > XiongZ_31763980.Prn-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-AlR_BETA,Prn-AlR_SE,Prn-AlR_P | sed -e 's/Prn-AlR_//g' | gzip -c > XiongZ_31763980.Prn-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-Sn_BETA,AlL-Sn_SE,AlL-Sn_P | sed -e 's/AlL-Sn_//g' | gzip -c > XiongZ_31763980.AlL-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-AlR_BETA,AlL-AlR_SE,AlL-AlR_P | sed -e 's/AlL-AlR_//g' | gzip -c > XiongZ_31763980.AlL-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Sn-AlR_BETA,Sn-AlR_SE,Sn-AlR_P | sed -e 's/Sn-AlR_//g' | gzip -c > XiongZ_31763980.Sn-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Ls-ChR_BETA,Ls-ChR_SE,Ls-ChR_P | sed -e 's/Ls-ChR_//g' | gzip -c > XiongZ_31763980.Ls-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Ls-Li_BETA,Ls-Li_SE,Ls-Li_P | sed -e 's/Ls-Li_//g' | gzip -c > XiongZ_31763980.Ls-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Ls-ChL_BETA,Ls-ChL_SE,Ls-ChL_P | sed -e 's/Ls-ChL_//g' | gzip -c > XiongZ_31763980.Ls-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ChR-Li_BETA,ChR-Li_SE,ChR-Li_P | sed -e 's/ChR-Li_//g' | gzip -c > XiongZ_31763980.ChR-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ChR-ChL_BETA,ChR-ChL_SE,ChR-ChL_P | sed -e 's/ChR-ChL_//g' | gzip -c > XiongZ_31763980.ChR-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Li-ChL_BETA,Li-ChL_SE,Li-ChL_P | sed -e 's/Li-ChL_//g' | gzip -c > XiongZ_31763980.Li-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-N_BETA,ExR-N_SE,ExR-N_P | sed -e 's/ExR-N_//g' | gzip -c > XiongZ_31763980.ExR-N.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-Prn_BETA,ExR-Prn_SE,ExR-Prn_P | sed -e 's/ExR-Prn_//g' | gzip -c > XiongZ_31763980.ExR-Prn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-AlL_BETA,ExR-AlL_SE,ExR-AlL_P | sed -e 's/ExR-AlL_//g' | gzip -c > XiongZ_31763980.ExR-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-Sn_BETA,ExR-Sn_SE,ExR-Sn_P | sed -e 's/ExR-Sn_//g' | gzip -c > XiongZ_31763980.ExR-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-AlR_BETA,ExR-AlR_SE,ExR-AlR_P | sed -e 's/ExR-AlR_//g' | gzip -c > XiongZ_31763980.ExR-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-N_BETA,EnR-N_SE,EnR-N_P | sed -e 's/EnR-N_//g' | gzip -c > XiongZ_31763980.EnR-N.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-Prn_BETA,EnR-Prn_SE,EnR-Prn_P | sed -e 's/EnR-Prn_//g' | gzip -c > XiongZ_31763980.EnR-Prn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-AlL_BETA,EnR-AlL_SE,EnR-AlL_P | sed -e 's/EnR-AlL_//g' | gzip -c > XiongZ_31763980.EnR-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-Sn_BETA,EnR-Sn_SE,EnR-Sn_P | sed -e 's/EnR-Sn_//g' | gzip -c > XiongZ_31763980.EnR-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-AlR_BETA,EnR-AlR_SE,EnR-AlR_P | sed -e 's/EnR-AlR_//g' | gzip -c > XiongZ_31763980.EnR-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-EnL_BETA,N-EnL_SE,N-EnL_P | sed -e 's/N-EnL_//g' | gzip -c > XiongZ_31763980.N-EnL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-EnL_BETA,Prn-EnL_SE,Prn-EnL_P | sed -e 's/Prn-EnL_//g' | gzip -c > XiongZ_31763980.Prn-EnL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-AlL_BETA,EnL-AlL_SE,EnL-AlL_P | sed -e 's/EnL-AlL_//g' | gzip -c > XiongZ_31763980.EnL-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-Sn_BETA,EnL-Sn_SE,EnL-Sn_P | sed -e 's/EnL-Sn_//g' | gzip -c > XiongZ_31763980.EnL-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-AlR_BETA,EnL-AlR_SE,EnL-AlR_P | sed -e 's/EnL-AlR_//g' | gzip -c > XiongZ_31763980.EnL-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-ExL_BETA,N-ExL_SE,N-ExL_P | sed -e 's/N-ExL_//g' | gzip -c > XiongZ_31763980.N-ExL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-ExL_BETA,Prn-ExL_SE,Prn-ExL_P | sed -e 's/Prn-ExL_//g' | gzip -c > XiongZ_31763980.Prn-ExL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-AlL_BETA,ExL-AlL_SE,ExL-AlL_P | sed -e 's/ExL-AlL_//g' | gzip -c > XiongZ_31763980.ExL-AlL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-Sn_BETA,ExL-Sn_SE,ExL-Sn_P | sed -e 's/ExL-Sn_//g' | gzip -c > XiongZ_31763980.ExL-Sn.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-AlR_BETA,ExL-AlR_SE,ExL-AlR_P | sed -e 's/ExL-AlR_//g' | gzip -c > XiongZ_31763980.ExL-AlR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-Ls_BETA,ExR-Ls_SE,ExR-Ls_P | sed -e 's/ExR-Ls_//g' | gzip -c > XiongZ_31763980.ExR-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-ChR_BETA,ExR-ChR_SE,ExR-ChR_P | sed -e 's/ExR-ChR_//g' | gzip -c > XiongZ_31763980.ExR-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-Li_BETA,ExR-Li_SE,ExR-Li_P | sed -e 's/ExR-Li_//g' | gzip -c > XiongZ_31763980.ExR-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExR-ChL_BETA,ExR-ChL_SE,ExR-ChL_P | sed -e 's/ExR-ChL_//g' | gzip -c > XiongZ_31763980.ExR-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-Ls_BETA,EnR-Ls_SE,EnR-Ls_P | sed -e 's/EnR-Ls_//g' | gzip -c > XiongZ_31763980.EnR-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-ChR_BETA,EnR-ChR_SE,EnR-ChR_P | sed -e 's/EnR-ChR_//g' | gzip -c > XiongZ_31763980.EnR-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-Li_BETA,EnR-Li_SE,EnR-Li_P | sed -e 's/EnR-Li_//g' | gzip -c > XiongZ_31763980.EnR-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnR-ChL_BETA,EnR-ChL_SE,EnR-ChL_P | sed -e 's/EnR-ChL_//g' | gzip -c > XiongZ_31763980.EnR-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-Ls_BETA,EnL-Ls_SE,EnL-Ls_P | sed -e 's/EnL-Ls_//g' | gzip -c > XiongZ_31763980.EnL-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-ChR_BETA,EnL-ChR_SE,EnL-ChR_P | sed -e 's/EnL-ChR_//g' | gzip -c > XiongZ_31763980.EnL-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-Li_BETA,EnL-Li_SE,EnL-Li_P | sed -e 's/EnL-Li_//g' | gzip -c > XiongZ_31763980.EnL-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,EnL-ChL_BETA,EnL-ChL_SE,EnL-ChL_P | sed -e 's/EnL-ChL_//g' | gzip -c > XiongZ_31763980.EnL-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-Ls_BETA,ExL-Ls_SE,ExL-Ls_P | sed -e 's/ExL-Ls_//g' | gzip -c > XiongZ_31763980.ExL-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-ChR_BETA,ExL-ChR_SE,ExL-ChR_P | sed -e 's/ExL-ChR_//g' | gzip -c > XiongZ_31763980.ExL-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-Li_BETA,ExL-Li_SE,ExL-Li_P | sed -e 's/ExL-Li_//g' | gzip -c > XiongZ_31763980.ExL-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,ExL-ChL_BETA,ExL-ChL_SE,ExL-ChL_P | sed -e 's/ExL-ChL_//g' | gzip -c > XiongZ_31763980.ExL-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-Ls_BETA,N-Ls_SE,N-Ls_P | sed -e 's/N-Ls_//g' | gzip -c > XiongZ_31763980.N-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-ChR_BETA,N-ChR_SE,N-ChR_P | sed -e 's/N-ChR_//g' | gzip -c > XiongZ_31763980.N-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-Li_BETA,N-Li_SE,N-Li_P | sed -e 's/N-Li_//g' | gzip -c > XiongZ_31763980.N-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,N-ChL_BETA,N-ChL_SE,N-ChL_P | sed -e 's/N-ChL_//g' | gzip -c > XiongZ_31763980.N-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-Ls_BETA,Prn-Ls_SE,Prn-Ls_P | sed -e 's/Prn-Ls_//g' | gzip -c > XiongZ_31763980.Prn-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-ChR_BETA,Prn-ChR_SE,Prn-ChR_P | sed -e 's/Prn-ChR_//g' | gzip -c > XiongZ_31763980.Prn-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-Li_BETA,Prn-Li_SE,Prn-Li_P | sed -e 's/Prn-Li_//g' | gzip -c > XiongZ_31763980.Prn-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Prn-ChL_BETA,Prn-ChL_SE,Prn-ChL_P | sed -e 's/Prn-ChL_//g' | gzip -c > XiongZ_31763980.Prn-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-Ls_BETA,AlL-Ls_SE,AlL-Ls_P | sed -e 's/AlL-Ls_//g' | gzip -c > XiongZ_31763980.AlL-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-ChR_BETA,AlL-ChR_SE,AlL-ChR_P | sed -e 's/AlL-ChR_//g' | gzip -c > XiongZ_31763980.AlL-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-Li_BETA,AlL-Li_SE,AlL-Li_P | sed -e 's/AlL-Li_//g' | gzip -c > XiongZ_31763980.AlL-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlL-ChL_BETA,AlL-ChL_SE,AlL-ChL_P | sed -e 's/AlL-ChL_//g' | gzip -c > XiongZ_31763980.AlL-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Sn-Ls_BETA,Sn-Ls_SE,Sn-Ls_P | sed -e 's/Sn-Ls_//g' | gzip -c > XiongZ_31763980.Sn-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Sn-ChR_BETA,Sn-ChR_SE,Sn-ChR_P | sed -e 's/Sn-ChR_//g' | gzip -c > XiongZ_31763980.Sn-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Sn-Li_BETA,Sn-Li_SE,Sn-Li_P | sed -e 's/Sn-Li_//g' | gzip -c > XiongZ_31763980.Sn-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,Sn-ChL_BETA,Sn-ChL_SE,Sn-ChL_P | sed -e 's/Sn-ChL_//g' | gzip -c > XiongZ_31763980.Sn-ChL.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlR-Ls_BETA,AlR-Ls_SE,AlR-Ls_P | sed -e 's/AlR-Ls_//g' | gzip -c > XiongZ_31763980.AlR-Ls.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlR-ChR_BETA,AlR-ChR_SE,AlR-ChR_P | sed -e 's/AlR-ChR_//g' | gzip -c > XiongZ_31763980.AlR-ChR.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlR-Li_BETA,AlR-Li_SE,AlR-Li_P | sed -e 's/AlR-Li_//g' | gzip -c > XiongZ_31763980.AlR-Li.csv.gz
zcat < XiongZ_31763980.txt.gz| csvcut -d " " -Sc CHR,BP,SNP,EA,OA,AlR-ChL_BETA,AlR-ChL_SE,AlR-ChL_P | sed -e 's/AlR-ChL_//g' | gzip -c > XiongZ_31763980.AlR-ChL.csv.gz
