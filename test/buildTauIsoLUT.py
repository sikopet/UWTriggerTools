import math
i = 0
isoflag=0
iso=0.2
offTauIso = 192.0
print"########################################\n# tauIsolation LUT for ISOL(A)= " + str(iso) +"  ISOL(B)= 100.00\n# Switch to ISOLB value at pt= "+str(offTauIso)+"\n# Address in binary: JJJJJJJJTTTTTTTT\n#<header> V1 16 1 </header>\n# Format:  Address  Payload  ## hwTauPt hwJetPt\n########################################"
while (i<65536):
	jetpt = math.floor(i/256)
	taupt = i-jetpt*256
	if (jetpt == 0 or taupt==0):
		jetpt = -1
		taupt = -1
		isoflag = 0
	else: 
		if ((jetpt-taupt)/taupt < iso or taupt >=offTauIso or jetpt >=255):
			isoflag = 1
		else:
			isoflag = 0

	print str(i) + " " + str(isoflag) + "\t## physPt tau: " + str(taupt) +" jet: "+str(jetpt)
	i = i+1



