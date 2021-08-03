#!/usr/bin/env python
# coding: utf-8





#--------------ΣΗΜΕΙΩΣΗ: ΤΑ FIGURES 1-27 ΑΝΤΙΣΤΟΙΧΟΥΝ ΣΤΟΝ ΑΜ el18175 ΕΝΩ TA 28-54 ΑΝΤΙΣΤΟΙΧΟΥΝ ΣΤΟΝ ΑΜ el18079





from numpy import pi,gcd,cos
from scipy import special
from sympy.combinatorics.graycode import GrayCode
import binascii
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.io.wavfile
import scipy.signal
import scipy.special


for program_counter in range(0,2):

		 #program_counter = 0 --> el18175
		 #program_counter = 1 --> el18079









	#---------------------------ΑΣΚΗΣΗ 1--------------------------------------










	AM = 5 if program_counter == 0 else 9;
	fm = 4000 if program_counter == 0 else 7000;

	F = np.gcd((AM+3)*fm,(AM+1)*fm); # βρίσκουμε την θεμελιώδη συχνότητα
	T = 1/F;

	# δειγματοληψία για fs1 = 20*fm

	fs1 = 20;   # σχετική συχνότητα δειγματοληψίας fs1/fm = 20
	t1 = np.linspace(0,4*T,4*fs1);      # 4*T για να έχουμε 4 περιόδους
	y1 = cos(2*pi*fm*t1)*cos(2*pi*(AM+2)*fm*t1);
	plt.figure();
	plt.plot(t1,y1,'.');
	plt.xlabel("t (sec)");
	plt.ylabel("Πλάτος (V)");
	plt.title("Δειγματοληπτημένο y(t) με fs1 = 20*fm");
	plt.grid(axis='y');

	# δειγματοληψία για fs2 = 100*fm
	fs2 = 100; # σχετική συχνότητα δειγματοληψίας fs1/fm = 100
	t2 = np.linspace(0,4*T,4*100);
	y2 = cos(2*pi*fm*t2)*cos(2*pi*(AM+2)*fm*t2);
	plt.figure();
	plt.plot(t2,y2,'.');
	plt.xlabel("t (sec)");
	plt.ylabel("Πλάτος (V)");
	plt.title("Δειγματοληπτημένο y(t) με fs1 = 100*fm");
	plt.grid(axis="y");

	# κοινή αναπαράσταση των παραπάνω δύο σημάτων
	plt.figure();
	plt.plot(t1,y1,'.',label="fs1 = 20 fm",color = "red");
	plt.plot(t2,y2,'.',label="fs2 = 100 fm",color = "green");
	plt.xlabel("t (sec)");
	plt.ylabel("Πλάτος (V)");
	plt.grid(axis="y");
	plt.title("Κοινό διάγραμμα για fs1 = 20 fm και fs2 = 100 fm");
	plt.legend(loc="upper right");

	# δειγματοληψία για fs = 5*fm
	fs = 5;   # σχετική συχνότητα δειγματοληψίας fs/fm = 5
	t3 = np.linspace(0,4*T,4*fs);
	y3 = cos(2*pi*fm*t3)*cos(2*pi*(AM+2)*fm*t3);
	plt.figure()
	plt.plot(t3,y3,'.');
	plt.xlabel("t (sec)");
	plt.ylabel("Πλάτος (V)");
	plt.grid();
	plt.title("Δειγματοληπτημένο y(t) με fs1 = 5 fm");












	# ----------------------------ΑΣΚΗΣΗ 2--------------------------------










	#παράμετροι
	fm=4000 if program_counter == 0 else 7000; #Hz
	AM=5 if program_counter == 0 else 9;
	A = 1   #V
	if (fm / 1000) % 2 == 1:
			bits = 5
	else:
			bits = 4
	Q = 2 * A / 2 ** bits + 1e-14

	def Mid_Rise_quantizer(x, Q):
					exit_Q = Q * (np.floor(y / Q) + .5)  #έξοδος κβαντιστή
					return exit_Q

	#a'
	def binary_to_gray(bits):
					bits = int(bits, 2)  # μετατροπή σε ακέραιο
					bits ^= (bits >> 1)
					return bin(bits)[2:].zfill(bits) #επιστροφή gray code


	#fm περιττή (βάσει του AM)
	if (bits == 5):
			k = '{0:05b}'.format(0)
			for x in range(1, 32, 1):
					k = np.append(k, '{0:05b}'.format(x))   #πίνακας k δυαδικών αριθμών της εξόδου του κβαντιστή exit_Q
			d = '{0:05b}'.format(0)
			for x in range(1, 32, 1):
					d = np.append(d, binary_to_gray(k[x]))  # πίνακας d με τον κωδικό Gray των k στοιχείων
			t1 = np.linspace(0, 4/fm, 4*20)  #fs=20fm
			y = np.cos(2 * np.pi * fm * t1) * np.cos(2 * np.pi * (2 + AM) * fm * t1) # είσοδος κβαντιστή
			exit_Q = Mid_Rise_quantizer(x, Q)  #κλήση του κβαντιστή

			#ορισμός άξονα y
			plt.figure();
			plt.yticks(np.arange(min(exit_Q), max(exit_Q) + 1e-14, Q),
										 ('00000', '00001', '00011', '00010', '00110', '00111', '00101', '00100', '01100',
											'01101', '01111', '01110', '01010', '01011', '01001', '01000', '11000', '11001',
											'11011', '11010', '11110', '11111', '11101', '11100', '10100', '10101', '10111',
											'10110', '10010', '10011', '10001', '10000'))
			plt.step(2 * t1, exit_Q, color='b', label='Κβαντισμένο σήμα')
			plt.axis()
			plt.xlabel('Χρόνος(sec)')
			plt.ylabel('Επίπεδα κβαντισμού')
			plt.title('Mid-Riser κβάντιση του σήματος y(t)')
			#plt.grid()
			plt.legend(loc='upper right')

	#fm άρτια (βάσει του ΑM)
	if (bits == 4):
			k = '{0:04b}'.format(0)
			for x in range(1, 16, 1):
					k = np.append(k, '{0:04b}'.format(x))  #πίνακας k δυαδικών αριθμών της εξόδου του κβαντιστή exit_Q
			d = '{0:04b}'.format(0)
			for x in range(1, 16, 1):
					d = np.append(d, binary_to_gray(k[x]))  # πίνακας d με τον κωδικό Gray των k στοιχείων
			t1 = np.linspace(0, 4/fm,4*20)  #fs=20fm
			y = np.cos(2 * np.pi * fm * t1) * np.cos(2 * np.pi * (2 + AM) * fm * t1) # είσοδος κβαντιστή
			exit_Q = Mid_Rise_quantizer(x, Q)  #κλήση του κβαντιστή

			#ορισμός άξονα y
			plt.figure();
			plt.yticks(np.arange(min(exit_Q), max(exit_Q) + 1e-14, Q),
										 ('0000', '0001', '0011', '0010', '0110', '0111', '0101', '0100',
											'1100', '1101', '1111', '1110', '1010', '1011', '1001', '1000'))
			plt.step(t1, exit_Q, color='b', label='Κβαντισμένο σήμα')
			plt.axis()
			plt.xlabel('Χρόνος(sec)')
			plt.ylabel('Επίπεδα κβαντισμού')
			plt.title('Mid-Riser κβάντιση του σήματος y(t)')
			#plt.grid()
			plt.legend(loc='upper right')



	#b'
	error=exit_Q-y #ορισμός σφάλματος
	#(i)
	std_10 = np.std(error[:10])  #Τυπική απόκλιση για τα πρώτα 10 δείγματα
	print('Τυπική απόκλιση για τα πρώτα 10 δείγματα: ', std_10)

	#(ii)
	std_20 = np.std(error[:20])  #Τυπική απόκλιση για τα πρώτα 20 δείγματα
	print('Τυπική απόκλιση για τα πρώτα 20 δείγματα: ', std_20)

	def dB(a):
			return 10*math.log10(a)

	#(iii)
	p1 = 0 #ισχύς
	p2 = 0
	np1 = 0 #ισχύς θορύβου
	np2 = 0

	# υπολογισμός ισχύων
	for i in range(10):
			p1 = p1 + y[i] ** 2
			p2 = p2 + y[i] ** 2
			np1 = np1 + error[i] ** 2
			np2 = np2 + error[i] ** 2
	for i in range(10, 20):
			p2 = p2 + y[i] ** 2
			np2 = np2 + error[i] ** 2
	#μέση τιμή ισχύων
	p1 = p1 / 10
	p2 = p2 / 20
	np1 = np1 / 10
	np2 = np2 / 20
	#υπολογισμός πειραματικού SNR
	SNR10 = dB(p1 / np1)  # SNR για τα πρώτα 10 δείγματα
	SNR20 = dB(p2 / np2)  # SNR για τα πρώτα 20 δείγματα
	print("SNR για τα πρώτα 10 δείγματα:", SNR10, "dB")
	print("SNR για τα πρώτα 20 δείγματα:", SNR20, "dB")
	#υπολογισμός θεωρητικού SNR
	SNR = dB((3/4)*2 ** (2 * bits))
	print("Θεωρητικό SNR:", SNR, "dB")

	#γ'
	# πίνακας τιμών Gray κωδικοποίησης συναρτήσει των bits
	if (bits == 5):
			grays = ['00000', '00001', '00011', '00010', '00110', '00111', '00101', '00100', '01100',
							 '01101', '01111', '01110', '01010', '01011', '01001', '01000', '11000', '11001',
							 '11011', '11010', '11110', '11111', '11101', '11100', '10100', '10101', '10111',
							 '10110', '10010', '10011', '10001', '10000']

	elif (bits == 4):
			grays = ['0000', '0001', '0011', '0010', '0110', '0111', '0101', '0100',
							 '1100', '1101', '1111', '1110', '1010', '1011', '1001', '1000']

	g = []
	x = np.zeros(20)
	levels = (np.arange(min(exit_Q), max(exit_Q) + 1e-14, Q))
	# εύρεση επιπέδων κβάντισης σε κωδικοποίηση Gray, για την έξοδο του κβαντιστή
	for i in range(20):
			x[i] = np.argmin(np.abs(levels - exit_Q[i]))

	for i in range(20):
			g.append(grays[int(x[i])])

	out = np.zeros(40*bits)
	bitstream = "".join(g)  # strings join --> πίνακας με bits
	# Polar RZ υλοποίηση: Πλάτη(0-->0, 1-->+A1,-A1 εναλλάξ)
	A1=fm/1000
	prev = -1
	for i in range(20*bits):
			if bitstream[i:(i+1)] == '1' and prev == -1:
					out[2*i] = A1
					prev = 1
			elif bitstream[i:(i+1)] == '1' and prev == 1:
					out[2*i] = -A1
					prev = -1
			else:
					out[2*i] = 0

	time = np.zeros(40*bits)
	for i in range(40*bits):
			time[i] = 0.0005*i  # Polar RZ υλοποίησηση: χρόνος
	plt.figure()
	plt.ylabel("Πλάτος[V]")
	plt.xlabel("Χρόνος[sec]")
	plt.title("Bitstream για μία περίοδο της εξόδου του κβαντιστή")
	plt.step(time, out, label='Volts που αναπαριστούν τα bits')
	plt.yticks([-A1, 0, A1])
	plt.grid()
	plt.legend(loc='upper right')













	#-----------------------------ΑΣΚΗΣΗ 3--------------------------------------














	length = 46;   # μήκος της ακολουθίας των bits
	Tb = 0.5;   # διάρκεια ενός bit

	A = 4 if program_counter == 0 else 7;

	seq = np.random.randint(0,2,length);  # η τυχαία ακολουθία bits
	t = np.linspace(0,length*Tb,length+1);  # ο άξονας του χρόνου
	t = np.append(0,t);   # προσθέτουμε ένα 0 στο διάνυσμα του χρόνου μόνο για λόγους σχεδίασης της γραφικής
	amp = [A if seq[i] == 1 else -A for i in range(0,length)];   # το διαμορφωμένο κατά B-PAM σήμα
	amp = np.append(0,amp);    # προσθέτουμε 0 στην αρχή
	amp = np.append(amp,0);    # και στο τέλος της κυματομορφής πάλι για λόγους σχεδίασης
	plt.figure();
	plt.step(t,amp,where = "post");
	plt.xticks(t[1::2]);
	plt.yticks(range(-A+1,A+1,1));
	plt.axhline(y=0, color='k', linewidth=0.7);
	plt.title("Τυχαία ακολουθία 46 bits διαμορφωμένη κατά B-PAM");
	plt.ylabel("Πλάτος (V)");
	plt.xlabel("t (sec)");
	plt.grid();

	Eb = (A**2)*Tb;     # μέση ενέργεια του bit
	s0 = -np.sqrt(Eb);  # αναπαράσταση της κυματομορφής που αντιστοιχεί στο 0 στο διάγραμμα αστερισμού
	s1 = np.sqrt(Eb);   # αναπαράσταση της κυματομορφής που αντιστοιχεί στο 1 στο διάγραμμα αστερισμού

	plt.figure();
	plt.plot([s0,s1],[0,0],'o');
	plt.axhline(y=0,linewidth=0.5,color='black');
	plt.axvline(x=0,linewidth=0.5,color='black');
	plt.annotate(0,[s0-0.1,0.005],size=12);
	plt.annotate(1,[s1,0.005],size=12);
	plt.title("Διάγραμμα αστερισμού τυχαίας ακολουθίας 46 bit διαμορφωμένης κατά B-PAM");
	plt.xlabel("I")
	plt.ylabel('Q')
	plt.axvline(x=0,color="black");
	plt.grid();

	reduced_amp = amp[1:-1];   # αφαιρούμε τα μηδενικά στην αρχή και στο τέλος γιατί δεν μας χρειάζονται πλέον

	t = np.arange(0,(length+1)*Tb,Tb);
	def increase_data_points(signal,factor):
			# αυτή η συνάρτηση επεκτείνει τον αριθμό των data points ενός υπάρχοντος σήματος κατά τον παράγοντα factor
			# για να μπορούμε να προσθέσουμε τον θόρυβο στην συνέχεια
			extended_signal = []
			for i in range(len(signal)):
					for j in range(factor):   # βάζουμε πολλές φορές την ίδια τιμή πχ αντί για [0,4] έχουμε [0,0,0,4,4,4]
							extended_signal.append(signal[i]);
			return extended_signal

	factor = 1000;  # παράγοντας αύξησης σημείων
	new_amp = increase_data_points(reduced_amp,factor);
	new_amp = np.concatenate([[0],new_amp,[0]]);
	new_t = np.linspace(Tb/factor,length*Tb,len(new_amp)-2);
	new_t = np.concatenate([[0],new_t,[23]]);

	#Eb/N0 = 5
	SNR = 5; # SNR σε dB
	N0 = Eb/(10**(SNR/10));

	#noise
	X_5dB = np.random.normal(0,np.sqrt(N0/2),len(new_amp)); # noise real part
	Y_5dB = np.random.normal(0,np.sqrt(N0/2),len(new_amp)); # noise imaginary part

	amp_noise_5dB = new_amp+X_5dB;   # το σήμα μετά την προσθήκη του πραγματικού μέρους του θορύβου
	plt.figure();
	plt.plot(new_t,amp_noise_5dB);
	plt.title("Διαμορφωμένη κατά BPAM κυματομορφή με προσθήκη AWGN, SNR = 5dB");
	plt.ylabel("Πλάτος (V)");
	plt.xlabel("t (sec)");
	plt.xticks(new_t[0::2*factor]);
	plt.yticks(range(-2*A,2*A,1 if A == 4 else 2));
	plt.grid();

	#Eb/N0 = 15
	SNR = 15; # SNR σε dB
	N0 = Eb/(10**(SNR/10));

	#noise
	X_15dB = np.random.normal(0,np.sqrt(N0/2),len(new_amp)); # noise real part
	Y_15dB = np.random.normal(0,np.sqrt(N0/2),len(new_amp)); # noise imaginary part

	amp_noise_15dB = new_amp+X_15dB;   # το σήμα μετά την προσθήκη του πραγματικού μέρους του θορύβου
	plt.figure();
	plt.plot(new_t,amp_noise_15dB);
	plt.title("Διαμορφωμένη κατά BPAM κυματομορφή με προσθήκη AWGN, SNR = 15dB");
	plt.ylabel("Πλάτος (V)");
	plt.xlabel("t (sec)");
	plt.xticks(new_t[0::2*factor]);
	plt.yticks(range(-2*A,2*A,1 if A == 4 else 2));
	plt.grid();

	Real = new_amp[factor::factor]/A*np.sqrt(Eb) + X_5dB[factor::factor];    # το πραγματικό μέρος στον αστερισμό
	# γράφουμε factor::factor για λόγους που έχουν να κάνουν με το πώς φτιάξαμε τον πίνακα των πλατών της κυματομορφής
	# 1. για να μην πάρουμε την αρχική τιμή που είναι 0 και υπάρχει μόνο για τον σκοπό της σχεδίασης
	# 2. για να πάρουμε ένα δείγμα ανά bit

	Im = 0+Y_5dB[factor::factor];    # το φανταστικό μέρος, όμοια με το πραγματικό

	plt.figure();
	plt.plot(Real,Im,'.',color='r');
	plt.plot([s0,s1],[0,0],'o',color='b');
	plt.axvline(x=0,color="black");
	plt.axis('equal');
	plt.xlim([-max(max(Real),abs(min(Real)))-0.1,max(max(Real),abs(min(Real)))+0.1]);  # για να είναι συμμετρικό ως προς τον κατακόρυφο άξονα
	plt.grid();
	plt.xlabel("I");
	plt.ylabel("Q");
	plt.title("Διάγραμμα αστερισμού B-PAM με προσθήκη AWGN, SNR = 5 dB");

	Real = new_amp[factor::factor]*np.sqrt(Eb)/A + X_15dB[factor::factor];
	# το πραγματικό μέρος στον αστερισμό όπως αναλύσαμε παραπάνω
	Im = 0+Y_15dB[factor::factor];    # το φανταστικό μέρος

	plt.figure();
	plt.plot(Real,Im,'.',color='r');
	plt.plot([s0,s1],[0,0],'o',color='b');
	plt.axvline(x=0,color="black");
	plt.axis('equal');
	plt.xlim([-max(max(Real),abs(min(Real)))-0.1,max(max(Real),abs(min(Real)))+0.1]);  # για να είναι συμμετρικό ως προς τον κατακόρυφο άξονα
	plt.grid();
	plt.xlabel("I");
	plt.ylabel("Q");
	plt.title("Διάγραμμα αστερισμού B-PAM με προσθήκη AWGN, SNR = 15 dB");

	bep_array = [];   # εδώ αποθηκεύουμε τα πειραματικά BER για SNR=0,...,15 dB
	bep_theoretical = [];   # αντίστοιχα τα θεωρητικά
	for snr_value in range(0,16):
			SNR = snr_value;
			N0 = Eb/(10**(SNR/10));    # προσαρμόζουμε τον θόρυβο για το εκάστοτε SNR
			length = 100000;   # το μήκος της συμβολοσειράς (μεγάλο για να συμφωνούν τα θεωρητικά με τα πειραματικά αποτελέσματα)

			bitstream = np.random.randint(0,2,size = length);  # τυχαία ακολουθία bits
			bitstream = [np.sqrt(Eb) if bitstream[i] > 0 else -np.sqrt(Eb) for i in range(0,length)];
			# διαμορφώνουμε κατά BPAM

			noise = np.random.normal(0,np.sqrt(N0/2),length); # πραγματικό μέρος του θορύβου
			bitstream_with_noise = bitstream+noise;  # προσθέτουμε τον θόρυβο
			decided = [np.sqrt(Eb) if bitstream_with_noise[i] > 0 else -np.sqrt(Eb) for i in range(0,length)];
			# η παραπάνω ακολουθία παράγεται με βάση την ευθεία απόφασης (κατακόρυφος άξονας)
			error = 0;    # αριθμός των λαθών
			for j in range(1,length-1):
					if decided[j] != bitstream[j]:   # αν διαφέρουν κάπου το εκπεμπόμενο με το ληφθέν σήμα, τότε έχει γίνει λάθος
							error += 1;

			bep_array.append(error/length);
			bep_theoretical.append(0.5*scipy.special.erfc(np.sqrt(Eb/N0)));

	plt.figure();
	plt.plot(bep_array,'o',color="purple",label = "Πειραματικό BER");
	plt.plot(bep_theoretical,label="Θεωρητικό BER");
	plt.title("Γραφική παράσταση του BER για SNR = 0...15dB");
	plt.legend();
	plt.xticks(range(0,16));
	plt.ylabel("Πιθανότητα λάθους (καθαρός αριθμός)");
	plt.xlabel("SNR (dB)");
	plt.grid();

	def PSK_seq(a,N,M):  #Δημιουργία ακολουθίας συμβόλων για MPSK από bitstream
			symbol_seq=[]
			k=int(math.log(M,2))
			for j in range(N//k):
					symbol=str(a[k*j])
					for i in range(1,k):
							symbol+=str(a[k*j+i])
					symbol_seq.append(symbol)
			return symbol_seq













	#----------------------------------ΑΣΚΗΣΗ 4----------------------------------














	#παράμετροι
	AM=5 if program_counter == 0 else 9
	A=4 if program_counter == 0 else 7
	fm=4000 if program_counter == 0 else 7000 #Hz
	AM_sum = 13 if program_counter == 0 else 16 # 1+7+5 = 13,   0+7+9 = 16

	Tb=0.5
	N1=10**5
	Es=(A**2)*Tb
	SNR=np.arange(0,16,1)
	Eb=Es/2
	N=46 #αριθμός δειγμάτων
	a=seq;  # η τυχαία ακολουθία από 46 bits
	qpsk=PSK_seq(a,N,4) #δημιουργία ακολουθίας συμβόλων για QPSK

	#a'
	def pi4_QPSK_constellation_diagram_impl(seq,E): #υλοποίηση διαγράμματος αστερισμού
			theta=math.pi/4
			IQ_stream=[]
			for x in seq:
					if x=="00":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*theta))
					elif x=="01":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*3*theta))
					elif x=="11":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*5*theta))
					elif x=="10":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*7*theta))
			return IQ_stream

	qpsk_constellation=pi4_QPSK_constellation_diagram_impl(qpsk,Es)
	plt.figure()
	plt.grid()
	plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Σήμα", color="green")
	plt.axvline(color="black",label="Περιοχή πόφασης")
	plt.axhline(color="black")
	plt.legend(loc=1,fontsize="small")
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.ylim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A))
	plt.xlim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A))
	plt.title("Διάγραμμα αστερισμού με κωδικοποίηση (π/4)Gray")
	plt.annotate("00",(math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("01",((-1)*math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("11",((-1)*math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));
	plt.annotate("10",(math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));

	#β'
	def AWGN_noise(Eb,ratio,N):
			N0=Eb/(10**(ratio/10))
			return np.random.normal(0,np.sqrt(N0/2),size=N)+1j*np.random.normal(0,np.sqrt(N0/2),size=N)

	plt.figure()
	plt.grid()
	r3=qpsk_constellation+AWGN_noise(Es,5,N//2)
	plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Μεταδιδόμενο σήμα",color="green")
	plt.scatter(np.real(r3),np.imag(r3),label="Λαμβανόμενο σήμα",color="red")
	plt.axvline(color="black",label="Περιοχή Απόφασης ")
	plt.axhline(color="black")
	plt.legend(loc=8,fontsize="small")
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.ylim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A))
	plt.xlim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A))
	plt.title("Διάγραμμα αστερισμού με κωδικοποίηση (π/4)Gray και πρόσθετο θόρυβο AWGN (SNR=5dB)")
	plt.annotate("00",(math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("01",((-1)*math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("11",((-1)*math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));
	plt.annotate("10",(math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));

	plt.figure()
	plt.grid()
	r4=qpsk_constellation+AWGN_noise(Es,15,N//2)
	plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Μεταδιδόμενο σήμα", color="green")
	plt.scatter(np.real(r4),np.imag(r4),label="Λαμβανόμενο σήμα",color="red")
	plt.axvline(color="black",label="Περιοχή Απόφασης ")
	plt.axhline(color="black")
	plt.legend(loc=1,fontsize="small")
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.ylim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A));
	plt.xlim((-1)*math.sqrt(Es)-math.sqrt(A),math.sqrt(Es)+math.sqrt(A));
	plt.title("Διάγραμμα αστερισμού με κωδικοποίηση (π/4)Gray και πρόσθετο θόρυβο AWGN (SNR=15dB)");
	plt.annotate("00",(math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("01",((-1)*math.sqrt(2*Es)/2,math.sqrt(2*Es)/2));
	plt.annotate("11",((-1)*math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));
	plt.annotate("10",(math.sqrt(2*Es)/2,(-1)*math.sqrt(2*Es)/2));

	#γ'
	def received_pi4_qpsk_impl(seq): #Αποδιαμόρφωση (π/4)-QPSK σήματος
			received=[]
			for j in range(len(seq)):
					if np.real(seq[j])>=0 and np.imag(seq[j])>=0:
							received.append("00")
					elif np.real(seq[j])<0 and np.imag(seq[j])>0:
							received.append("01")
					elif np.real(seq[j])<=0 and np.imag(seq[j])<0:
							received.append("11")
					elif np.real(seq[j])>0 and np.imag(seq[j])<0:
							received.append("10")
			return received

	def reconcil_qpsk(seq,received): #Υπολογισμός του αριμθού και του ρυθμού σφαλμάτων
			count=0
			for j in range(len(seq)):
					if seq[j][0]!=received[j][0]:
							count+=1
					if seq[j][1]!=received[j][1]:
							count+=1
			return count, count/(2*len(seq))

	c=np.random.randint(2,size=N1)
	qpsk1=PSK_seq(c,N1,4)
	constellation_qpsk1=pi4_QPSK_constellation_diagram_impl(qpsk1,Es)
	BER2=[]
	BER2_th=[]

	for ratio in range(16):
			r5=constellation_qpsk1+AWGN_noise(Es,ratio,N1//2)
			received=received_pi4_qpsk_impl(r5)
			er,BER_=reconcil_qpsk(qpsk1,received)
			BER2.append(BER_)
			BER2_th.append(0.5*special.erfc(math.sqrt(0.5*10**(ratio/10))))

	print ("Πειραματική τιμή BER:",BER2)
	print ("Θεωρητική τιμή BER:",BER2_th)

	fig, (ax1,ax2) = plt.subplots(2,1,sharex=True, sharey=True)
	ax1.set_title("BER diagram")
	ax1.plot(SNR,BER2,"gx--",linewidth=0.5,label="Πειραματική τιμή BER")
	ax1.set_ylabel("BER")
	ax1.grid()
	ax1.legend()
	ax2.plot(SNR,BER2_th,"yx--",linewidth=0.5,label="Θεωρητική τιμή BER")
	plt.xlabel("SNR $E_s/N_0$ (dB)")
	plt.ylabel("BER")
	plt.grid()
	plt.legend()

	#δ', (i-iii)
	def text_to_bits(text, encoding='utf-8'):
			bits = bin(int.from_bytes(text.encode(encoding), 'big'))[2:]
			return bits.zfill(8 * ((len(bits) + 7) // 8))

	def text_from_bits(bits):
			q="".join(bits)
			bytes_=[q[8*i:8*i+8] for i in range(len(q)//8)] #διαχωρισμός του bitsream σε bytes(8bits)
			return "".join([chr(int(b,2)) for b in bytes_])

	filename="shannon_even.txt" if AM_sum%2==0 else "shannon_odd.txt"
	f=open(filename,"r")
	text=f.read()
	f.close()
	k=text_to_bits(text)
	k=k.replace("0b","")

	#(iv-v)
	qpsk_txt=PSK_seq(k,len(k),4)
	Es_txt=Tb
	Eb_txt=Es_txt/2

	def QPSK_constellation_diagram_impl(seq,E):
			IQ_stream=[]
			for x in seq:
					if x=="00":
							IQ_stream.append(math.sqrt(E)*np.exp(0))
					elif x=="01":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*math.pi/2))
					elif x=="11":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*math.pi))
					elif x=="10":
							IQ_stream.append(math.sqrt(E)*np.exp(1j*3*math.pi/2))
			return IQ_stream

	x=[-5,5]
	y1,y2=[-5,5],[5,-5]
	qpsk_txt_constellation=QPSK_constellation_diagram_impl(qpsk_txt,Es_txt)
	plt.figure()
	plt.annotate("00",(math.sqrt(Es_txt),0))
	plt.annotate("01",(0,math.sqrt(Es_txt)))
	plt.annotate("11",((-1)*math.sqrt(Es_txt),0))
	plt.annotate("10",(0,(-1)*math.sqrt(Es_txt)))
	plt.grid()
	plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation),color="green")
	plt.plot(x,y1,color="black",label="Περιοχή Απόφασης ")
	plt.plot(x,y2,color="black")
	plt.legend(loc=0,fontsize="small")
	plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.title("Διάγραμμα αστερισμού για QPSK διαμόρφωση του αρχείου κειμένου")

	r6=qpsk_txt_constellation+AWGN_noise(Es_txt,5,len(qpsk_txt))
	plt.figure()
	plt.annotate("00",(math.sqrt(Es_txt),0))
	plt.annotate("01",(0,math.sqrt(Es_txt)))
	plt.annotate("11",((-1)*math.sqrt(Es_txt),0))
	plt.annotate("10",(0,(-1)*math.sqrt(Es_txt)))
	plt.grid()
	plt.scatter(np.real(r6),np.imag(r6),label="Λαμβανόμενο σήμα",s=15,color="red")
	plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation),label="Μεταδιδόμενο σήμα",s=15,color="green")
	plt.plot(x,y1,color="black",label="Περιοχή Απόφασης ")
	plt.plot(x,y2,color="black")
	plt.legend(loc=0,fontsize="small")
	plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.title("Διάγραμμα αστερισμού για QPSK διαμόρφωση του αρχείου κειμένου με πρόσθετο θόρυβο AWGN 5dB")

	plt.figure()
	r7=qpsk_txt_constellation+AWGN_noise(Es_txt,15,len(qpsk_txt))
	plt.annotate("00",(math.sqrt(Es_txt),0))
	plt.annotate("01",(0,math.sqrt(Es_txt)))
	plt.annotate("11",((-1)*math.sqrt(Es_txt),0))
	plt.annotate("10",(0,(-1)*math.sqrt(Es_txt)))
	plt.grid()
	plt.scatter(np.real(r7),np.imag(r7),label="Λαμβανόμενο σήμα",s=15,color="red")
	plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation),label="Μεταδιδόμενο σήμα",s=15,color="green")
	plt.plot(x,y1,color="black",label="Περιοχή Απόφασης ")
	plt.plot(x,y2,color="black")
	plt.legend(loc=0,fontsize="small")
	plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
	plt.xlabel("I")
	plt.ylabel("Q")
	plt.title("Διάγραμμα αστερισμού για QPSK διαμόρφωση του αρχείου κειμένου με πρόσθετο θόρυβο AWGN 15dB");

	def received_qpsk_impl(seq):
			received=[]
			for j in range(len(seq)):
					if abs(np.real(seq[j]))>abs(np.imag(seq[j])):
							if np.real(seq[j])>0:
								 received.append("00")
							else:
								 received.append("11")
					else:
							if np.imag(seq[j])>0:
								 received.append("01")
							else:
								 received.append("10")
			return received

	#vi
	r6_=received_qpsk_impl(r6)
	be_r6,BER_r6=reconcil_qpsk(qpsk_txt,r6_)
	print("Πειραματική τιμή BER με SNR={}dB: {}({} σφάλματα σε {}bits)".format(5,BER_r6,be_r6,len(k)))
	print("Θεωρητική τιμή BER με SNR={}dB: {}".format(5,0.5*special.erfc(math.sqrt(0.5*10**(5/10)))))

	r7_=received_qpsk_impl(r7)
	be_r7,BER_r7=reconcil_qpsk(qpsk_txt,r7_)
	print("Πειραματική τιμή BER με SNR={}dB: {}({} σφάλματα σε {} bits)".format(15,BER_r7,be_r7,len(k)))
	print("Θεωρητική τιμή BER με SNR={}dB: {}".format(15,0.5*special.erfc(math.sqrt(0.5*10**(15/10)))))

	#vii
	def reconstruct_text(received,SNR):
			received_bitstream=[]
			for c in received:
					received_bitstream+=c
			received_bitstream="".join(received_bitstream)
			reconstructed_text=text_from_bits(received_bitstream)
			print("Ανακατασκευασμένο αρχείο κειμένου με {}dB θόρυβο AWGN:".format(SNR))
			print(reconstructed_text)
			f_w=open(filename.replace(".txt","")+' with noise {}dB.txt'.format(SNR),'w',encoding='utf-8')
			f_w.write(reconstructed_text)
			f_w.close()
	reconstruct_text(r6_,5)
	print(" ")
	reconstruct_text(r7_,15)














	#------------------------------ΑΣΚΗΣΗ 5-------------------------------------



















	soundfile = "soundfile1_lab2.wav" if AM_sum % 2 == 1 else "soundfile2_lab2.wav";

	rate,data = scipy.io.wavfile.read(soundfile); # διάβασμα του αρχείου
	t = np.linspace(0,len(data),len(data));
	t *= 1/rate;  # κανονικοποίηση του άξονα του χρόνου για να δείχνει second
	plt.figure();
	plt.plot(t,data);
	plt.xlabel("t (sec)");
	plt.title("Κυματομορφή του .wav αρχείου που διαβάστηκε");
	plt.grid();

	R = 8;   # bits κβάντισης
	L = 2**R;    # στάθμες κβάντισης

	max1 = max(max(data),abs(min(data)));
	D = 2*max1/L;   # μέγεθος βήματος

	def mid_riser(x):    #  παίρνει την είσοδο x και επιστρέφει την έξοδο του κβαντιστή
			return D*(math.floor(x/D)+1/2);

	i = mid_riser(min(data)) # χαμηλότερο επίπεδο κβάντισης
	quantization_levels = [];
	for j in range(0,L):     # δημιουργώ τα επίπεδα της κβάντισης
			quantization_levels.append(round(i,6))
			i += D

	quantized_signal = list(map(mid_riser,data))  # εκτελούμε την κβάντιση του σήματος
	for i in range(len(quantized_signal)):     # στρογγυλοποιήσεις για να έχουμε ακριβή αντιστοιχία με τον κβαντιστή
			quantized_signal[i] = round(quantized_signal[i],6)
	gray = GrayCode(R)
	gray = list(gray.generate_gray()); # κώδικας Gray για τα επίπεδα κβάντισης


	plt.figure();
	plt.yticks(quantization_levels[0::10]);
	plt.plot(t,quantized_signal,color="darkred");  # η έξοδος του κβαντιστή
	plt.grid();
	plt.title("Σήμα που προκύπτει από κβάντιση Mid Riser 8-bit");
	plt.ylabel("Έξοδος κβαντιστή");
	plt.xlabel("t (sec)");

	def gray_quantization(quantized_signal,gray,quantization_levels):  # αντιστοιχεί τιμές του κβαντισμένου σήματος σε κώδικα Gray
			array = [];
			for i in range(0,len(quantized_signal)):
					y = quantization_levels.index(quantized_signal[i]);
					array.append(gray[y]);
			return array;

	array =gray_quantization(quantized_signal,gray,quantization_levels);
	s = "".join(array);  # εν΄ώνουμε τις 8-΄άδες που προέκυψαν σε μία ενιαία ακολουθία bits


	def generateMPSK(s,length,M): # δημιουργεί ακολουθία συμβόλων M-PSK λαμβάνοντας ως είσοδο μία ακολουθία bits
			sequence = []
			k = int(math.log(M,2));
			for i in range(0,length//k):
					symbol = str(s[k*i]);
					for j in range(1,k):
							symbol += str(s[k*i+j]);
					sequence.append(symbol)
			return sequence

	qpsk_signal = generateMPSK(s,len(s),4);

	Tb = 0.5; # διάρκεια συμβόλου

	Es = Tb; # ενέργεια συμβόλου
	Eb = Es/2;
	def constellation_plotter():    # δημιουργεί τα σημεία 00 01 11 10 και τους βάζει ταμπέλες
			plt.annotate("00",(np.sqrt(Es),0));
			plt.annotate("01",(0,np.sqrt(Es)));
			plt.annotate("11",(-np.sqrt(Es),0));
			plt.annotate("10",(0,-np.sqrt(Es)));
			plt.axline([0,0],[1,1],color='black');
			plt.axline([0,0],[-1,1],color='black');
			plt.plot(np.sqrt(Es),0,'o',0,np.sqrt(Es),'o',-np.sqrt(Es),0,'o',0,-np.sqrt(Es),'o',color="red");
			plt.plot(np.sqrt(Es),0,'o',color="red",label="εκπεμπόμενο σήμα");

	def noise(E,SNR,length):   # δημιουργεί τον AWGN με κατάλληλο N0 κάθε φορά, ανάλογα το SNR
			N0 = E/(10**(SNR/10))
			noise = np.random.normal(0,np.sqrt(N0/2),size=length)+1j*np.random.normal(0,np.sqrt(N0/2),size=length)
			return noise

	plt.figure();
	constellation_plotter();
	plt.xlabel("I");
	plt.ylabel("Q");
	plt.title("Διάγραμμα αστερισμού σήματος QPSK");
	plt.xlim([-2*Es,2*Es]);
	plt.ylim([-2*Es,2*Es]);
	plt.grid();

	def constellation_from_data(data,E):
			# δημιουργεί την απεικόνιση του κάθε ληφθέντος συμβόλου στον αστερισμό
			# ώστε μετά να μπορούμε να προσθέσουμε τον θόρυβο στη συνέχεια
			s = []
			for b in data:
					if b == "00":
							s.append(np.sqrt(E)*np.exp(0));
					elif b == "01":
							s.append(np.sqrt(E)*np.exp(1j*pi/2))
					elif b == "11":
							s.append(np.sqrt(E)*np.exp(1j*pi))
					elif b == "10":
							s.append(np.sqrt(E)*np.exp(3j*pi/2))
			return s

	initial_signal = constellation_from_data(qpsk_signal,Es) # το αρχικό σήμα
	noisy_signal_4dB = initial_signal + noise(Es,4,len(qpsk_signal)); # το αρχικό σήμα μαζί με θόρυβο (SNR = 4dB)
	noisy_signal_14dB = initial_signal + noise(Es,14,len(qpsk_signal)); # το αρχικό σήμα μαζί με θόρυβο (SNR = 14dB)

	#εδώ απεικονίζουμε τα διαγράμματα αστερισμού μετά από προσθήκη θορύβου (SNR =4dB, 14dB αντίστοιχα)
	plt.figure();
	plt.plot(np.real(noisy_signal_4dB),np.imag(noisy_signal_4dB),'.',color = "green",label='Ληφθέν σήμα μετά από θόρυβο');
	constellation_plotter();
	plt.legend(loc = "upper right");
	plt.xlabel("I");
	plt.ylabel("Q");
	plt.grid();
	plt.title("Διάγραμμα αστερισμού QPSK μετά από προσθήκη θορύβου, SNR = 4dB")

	plt.figure();
	plt.plot(np.real(noisy_signal_14dB),np.imag(noisy_signal_14dB),'.',color = "green",label='Ληφθέν σήμα μετά από θόρυβο');
	constellation_plotter();
	plt.legend(loc = "upper right");
	plt.xlabel("I");
	plt.ylabel("Q");
	plt.grid();
	plt.title("Διάγραμμα αστερισμού QPSK μετά από προσθήκη θορύβου, SNR = 14dB")

	def QPSK_decoder(stream):  # εκτελεί την αποδιαμόρφωση του σήματος κατά QPSK
			decoded = []
			for i in range(0,len(stream)):
					x = np.real(stream[i]);   # πραγματικό μέρος
					y = np.imag(stream[i]);   # και φανταστικό μέρος, βλέπουμε σε ποια από τις 4 περιοχές απόφασης ανήκουν
					if x > np.abs(y): decoded.append("00");  # περιοχή 00
					elif y > np.abs(x): decoded.append("01");  # περιοχή 01
					elif x < -np.abs(y): decoded.append("11"); # περιοχή 11
					elif y < -np.abs(x): decoded.append("10");   # περιοχή 10
			return decoded

	decoded_4dB = QPSK_decoder(noisy_signal_4dB);
	decoded_14dB = QPSK_decoder(noisy_signal_14dB);

	def error_rate(original,received):    # πειραματικός υπολογισμός του BER
			errors = 0;
			for i in range(len(original)):
					if original[i][0] != received[i][0]:   # αν γίνει λάθος στο πρώτο bit του συμβόλου
							errors += 1;
					if original[i][1] != received[i][1]:   # αν γίνει λάθος στο δεύτερο bit του συμβόλου
							errors += 1;
			return errors/(2*len(original));     # διαιρούμε με το πλήθος των bits (2*αριθμός συμβόλων) γιατί θέλουμε πιθανότητα

	error_rate_4dB = error_rate(qpsk_signal,decoded_4dB);   # πειραματικό BER (SNR = 4dB)
	error_rate_14dB = error_rate(qpsk_signal,decoded_14dB);   # πειραματικό BER (SNR = 14dB)
	theoretical_BER_4dB = 0.5*special.erfc(np.sqrt(0.5*10**(4/10))); # θεωρητικό BER (SNR = 4dB)
	theoretical_BER_14dB = 0.5*special.erfc(np.sqrt(0.5*10**(14/10))); # θεωρητικό BER (SNR = 14dB)

	print("Πειραματικό BER με SNR=4dB : {}".format(error_rate_4dB));
	print("Θεωρητικό BER με SNR=4dB : {}".format(theoretical_BER_4dB))
	print("Πειραματικό BER με SNR=14dB : {}".format(error_rate_14dB));
	print("Θεωρητικό BER με SNR=14dB : {}".format(theoretical_BER_14dB))

	#γράφουμε το αρχείο και πλοτάρουμε το αρχείο που γράψαμε, για 4dB και 14dB αντίστοιχα
	def write_soundfile(received,SNR):
			received_=[];
			for i in range(0,len(received)//4):
					received_.append(received[4*i]+received[4*i+1]+received[4*i+2]+received[4*i+3]);  # φτιάχνουμε οκτάδες
			u = []
			for w in received_:
					u.append(gray.index(w));  # αντιστοιχίζουμε τις οκτάδες στις τιμές 0-255
			t = np.linspace(0,len(data),len(data));
			t *= 1/rate;
			u = np.array(u);
			u = u.astype("uint8")  # μετατρέπουμε τα δεδομένα σε 8-bit unsigned
			plt.figure();
			plt.plot(t,u);
			plt.xlabel("t (sec)");
			plt.title("Ανακατασκευασμένο σήμα μετά από προσθήκη θορύβου με SNR = {} dB".format(SNR));
			plt.ylabel("Κανονικοποιημένη έξοδος κβαντιστή");
			plt.grid();
			file_no_extension = "soundfile1_lab2_" if soundfile == "soundfile1_lab2.wav" else "soundfile2_lab2_"  # remove extension
			scipy.io.wavfile.write(file_no_extension+"with_noise_{}dB.wav".format(SNR),rate,u);

	write_soundfile(decoded_4dB,4);
	write_soundfile(decoded_14dB,14);

plt.show();