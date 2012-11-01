'''
Last update 20110501. This works with python version 2.7.2 or lower.
'''

def r2b(n,b):
	'''
	It will return the base b representation of real number, base 10, n.
	'''
	print 'First the whole part'
	r = int(n)
	f = n -r
	L = []
	D = []
	while (r != 0):
		x = r
		r, a =  r/b, r%b
		print x , ' = ', r, '*', b, ' + ', a
		L.insert(0,a)           #could have used expression for "a" directly
	s = ''
	print 'Now, the fraction part'
	for x in L:
		s = s + str(x)          #just to do a pretty print out          
	for i in range(20):
		d = b*f
		print d, ' = ', f, '*', b
		D.append(int(d))
		f = d - int(d)
	s = s + "."
	for x in D:
		s = s + str(x)
	print
	print s
	

def dc(n,b):
	'''
	It will convert decimal n to base b, 2 <=b <=9.
	'''
	L = []
	if b >=10: 
		return ("Can't handle b >= 10 for now.")
	while (n != 0):
		a = n%b  			#remainder
		n = n/b				# quotient
		L.insert(0,a)		#could have used expression for "a" directly
	s = ''
	for x in L:
		s = s + str(x)		#just to do a pretty print out
	print s
		
def sqrt(n,p):
	'''
	It will return the square root of n with precision p.
	'''
	A = 50
	B = n/A
	k = 0
	print "Root A is ", A, "\t", "Root B is ", B, "\t", "Iteration is ", k
	while (abs(A-B) >= p):
		A = (A+B)/2.
		B = n/A
		k = k+1
		print "Root A is ", A, "\t", "Root B is ", B, "\t", "Iteration is ", k
	return B

def nthroot(A, n, p):
        '''
        It returns the nth root of A with precision p
        '''
        x = 1.           #initial guess
        y = A
        while (abs(x-y) > p):
                x = (1./n)*((n-1)*x + y)
                y = A/(x**(n-1.))
        return y
def isA(x):
	'''
	Is this an algorithm?
	'''
	L=[]
	L.append(x)
	while (x !=1):
		if x%2 == 0:
			x = x/2
			L.append(x)
		else:
			x = 3*x + 1
			L.append(x)
	return L
	

def isbn2():
	'''
	It will print the ISBN checksum
	'''
	s = raw_input("Please enter a 9-digit ISBN ")
	i = 1
	x = 0
	for c in s:
		x = x + i * int(c)
		i = i+1
	cd = x%11
	if cd == 10:
		print "The check digit is x."
	else:
		print "The check digit is ", cd
	
def isbn(s = ''):
	'''
	It will print the ISBN checksum
	'''
	i = 1
	x = 0
	for c in s:
		x = x + i * int(c)
		i = i+1
	cd = x%11
	if cd == 10:
		print "The check digit is x."
	else:
		print "The check digit is ", cd
	
def hex(n):
	'''
	It will return the hexidecimal representation of n.
	'''
	L = []
	while(n !=0):
		a = n%16  				#remainder
		n = n/16				# quotient
		L.insert(0,a)			#could have used expression for "a" directly
	s = ''
	for x in L:
		if x == 10:
			c = 'A'
		if x == 11:
			c = 'B'
		if x == 12:
			c = 'C'
		if x == 13:
			c = 'D'
		if x == 14:
			c = 'E'
		if x == 15:
			c = 'F'
		if x < 10:
			c = str(x)
		s = s + c			
	return s	
			

def prime_factors(n):
	'''
	It will return the prime factorization of n.
	'''
	import math
	i = 2
	L = []
	s = str(n)+'=1'
	x = math.sqrt(n)
	while n!=1:
		if n <= 0: break 		#To prevent infinite loop by bad input
		if i > x:
			s = s + '*' + str(n)
			break
		if n%i == 0:
			#L.append(i)
			s = s+'*'+str(i)
			n = n/i
		else:
			i=i+1
	return s


def primes(n):
	'''
	It writes to a file all primes up to n.
	'''
	import sys
	pfile = open('C:\Users\Esfahanian\Desktop\CourseWare\CSE201_Spring_2010\Homework\Primes.txt','a')
	for i in range(n):
		s = prime_factors(i)
		if s.count('*') == 1:
			print >>pfile, i
	pfile.close

def nextprime(n=2):
	'''
	It outputs the next prime.
	'''
	while 1:
		s = prime_factors(n)
		if s.count('*') ==1:
			n = n + 1
			yield n - 1
		else:
			n = n + 1

def random_list(n, m):
	'''
	It will return a list of n random numbers between 0 and m.
	'''
	L = []
	import random
	for i in range(n):
		L.append(random.choice(range(m)))
	return L

def exchange_sort(L = [8, 24, 21, 2, 17, 22, 7]):
        print 'This is pass #', 0, L
        for i in range(len(L)-1):
                for j in range(i, len(L)):
                        if L[j] < L[i]:
                                a = L[j]
                                L[j] = L[i]
                                L[i] = a
                print 'This is pass #', i+1, L

def bubble_sort(L = [2, 7, 8, 17, 21, 22, 24]):
        print 'This is pass #', 0, L
        for i in range(len(L)-1):
                swaps = 0
                for j in range(len(L)-i-1):
                        if L[j] > L[j+1]:
                                a = L[j+1]
                                L[j+1] = L[j]
                                L[j] = a
                                swaps = swaps + 1
                print 'This is pass #', i+1, L
                if swaps == 0:
                        break

def quick_sort(L=[]):
	'''
	It will sort L using quicksort algorithm.
	'''
	if len(L) <= 1:
		return L
	else:
		import random
		import copy
		l = copy.deepcopy(L)	# To keep the input list intact
		print "The working list is ", l
		left = []
		right = []
		m = random.choice(range(len(l)))
		pivot = l[m]
		print "The pivot is ", pivot
		del l[m]
		for x in l:
			if x < pivot:
				left.append(x)
			else:
				right.append(x)
		print "The left list is ", left
		print "The right list is ", right
		print "-------- end of pass ---------"
		return quick_sort(left) + [pivot] + quick_sort(right)

#another quicksort		

def qsort(list):
    if len(list) <= 1:
        return list
    else:
        pivot = list[0]
        less, greater, equal = [], [], []
        for x in list:
            if x < pivot:
                less.append( x )
            elif x == pivot:
                equal.append( x )
            else:
                greater.append( x )

    return qsort(less) + equal + qsort(greater)
		
def easy_sort(L=[]):
	'''
	It will sort L using a simple recursive sort algorithm.
	'''
	if len(L) <= 1:
		return L
	else:
		import copy
		l = copy.deepcopy(L)		# To keep the input list intact
		min_value = min(l)
		m = l.index(min_value)
		del l[m]
		return [min_value] + easy_sort(l)
		
def Fourier_square_wave(n = 1):
	'''
	It will plot n harmonics of the square wave.
	'''
	import pylab
	import numpy
	t = numpy.arange(0,1,0.001)
	s = 0
	for i in range(1,n,2):
		s = s + (1/(i*numpy.pi))*numpy.sin(2*numpy.pi*i*t)
	pylab.plot(t,s,lw=3)
	pylab.show()

def Fourier_triangle_wave(n = 1):
	'''
	It will plot n harmonics of the triangle wave.
	'''
	import pylab
	import numpy
	t = numpy.arange(0,2,0.001)
	s = 0
	pi = numpy.pi
	for i in range(1,n,2):
		s = s + (8/(pi**2))*((-1)**((i-1)/2))*(numpy.sin(pi*i*t))/(i**2)
	pylab.grid(True)
	pylab.plot(t,s,lw=3)
	pylab.show()

def palindrome(n):
	'''
	It will return n (even) length palindromes.
	'''
	L = []
	for i in range(10**(n/2 -1), 10**(n/2), 1):
		s = str(i)
		b = ""
		for j in range(len(s)):
			b = b + s[len(s) - j - 1]
		p = int(s+b)
		L.append(p)
	return L
	
def bd(k,n =23):
	'''
	This is to verify the claim that if there are 23 people in a room there is more than 50% likelihood
	that two of them would have the same birthday where birthday is any of the 365 days.
	'''
	import random
	hit = 0
	e = 0.0
	for x in range(k):
		L = []
		for i in range(n):
			d = random.randint(1,365)
			L.append(d)
		#L.sort()
		#print " L = ", L
		if len(L) > len(set(L)):
			hit = hit + 1
		#print "# of Hits = ", hit
	print "hit prob = ", hit/(k+0.)
			
def mytest():
	import pylab
	import numpy
	pi = numpy.pi
	t = numpy.arange(0,4*pi,0.0001)
	g = numpy.sin(3*t) + numpy.sin(4*t)			# the blow curve
	h = numpy.sin(t)				# the green curve
	pylab.grid(True)
	pylab.plot(t,g,t,h, lw=3)
	pylab.show()
	
def periodic():
	import pylab
	import numpy
	pi = numpy.pi
	T = 2*pi/37
	t = numpy.arange(0,2*T,0.0001)
	g = numpy.sin(37*t)			# the blow curve
	h = numpy.sin(37*(t+T))							# the green curve
	pylab.grid(True)
	#pylab.hold(True)
	pylab.subplot(211)
	pylab.plot(t,g, 'g',  lw=3)
	pylab.subplot(212)
	pylab.plot(t,h, 'b',  lw=3)
	pylab.show()
	
def cc(s, shift):
	s= s.upper()
	import string
	L = string.ascii_uppercase
	code = ''
	for c in s:
		if c in L:
			i = L.index(c)
			j = (i + shift)%26
			code = code + L[j]
		else:
			code = code + c
	print code.lower()

def vig_enc(p = 'HELLO', k = 'XMCKL'):
	'''
	p is the mesage to be encrypted, k is the key. 
	'''
	
	import string
	L = string.ascii_uppercase

	p = p.upper()
	k= k.upper()
	code = ''
	from math import ceil
	r = int(ceil(len(p)/len(k)))
	#print r
	for i in range(r+1):
		k = k + k
	k = k[0:len(p)]
	print p
	print k
	for i in range(len(p)):
		c = p[i]
		r = k[i]
		i = L.index(c)
		j = L.index(r)
		j = (i + j)%26
		code =  code + L[j]
	print code


def vig_dec(p = 'EQNVZ' , k = 'TQURI'):

	'''
	p is the encrypted message, and k is the key
	'''       
	import string
	L = string.ascii_uppercase

	p = p.upper()
	k= k.upper()
	code = ''
	from math import ceil
	r = int(ceil(len(p)/len(p)))
	#print r
	for i in range(r+1):
		k = k + k
	k = k[0:len(p)]
	print p
	print k
	for i in range(len(k)):
		c = k[i]
		r = p[i]
		i = L.index(c)
		j = L.index(r)
		j = (j - i)%26
		code =  code + L[j]
	print code

def BA(n,m,b):
	'''
	Adds n and m in base b and returns the result
	'''
	if n >= m:
		n1 = n
		n2 = m
	else:
		n1 = m
		n2 = n
	L1 = []
	L2 = []
	sum = []
	for c in str(n1):
		L1.insert(0,int(c))
	for c in str(n2):
		L2.insert(0,int(c))
	s = ''
	for i in range(len(L1) - len(L2)):
		L2.append(0)
	#print L1, L2
	carry = 0
	for x in range(len(L1)):
		d = (L1[x] + L2[x] + carry)%b
		carry = (L1[x] + L2[x] + carry)/b
		s = str(d) + s
		sum.insert(0,d)
	if carry != 0:
		s = str(carry) + s
		sum.insert(0,carry)
	#print int(s)
	L1.reverse()
	L2.reverse()
	print ' ',
	for x in L1:
		print x,
	print
	print '+',
	for x in L2:
		print x,
	print
	line = '-'
	line = line *(2*len(L1) + 1)
	print line
	if len(sum) > len(L1):
		for x in sum:
			print x,
	else:
		print ' ',
		for x in sum:
			print x,		
	print
	
def CharCount(s =''):
	'''
	Finds the frequency of each character in string s
	'''
	l = []
	k = []
	s = s.lower()
	for char in s:
		if char != " " and char != "." and char not in l:
			l.append(char)
	l.sort()
	for char in l:
		k.append((s.count(char), char))
		print "The count of", char, " is" , char*s.count(char)
	print "The most frequent character is **", max(k)[1], "** with frequency ", max(k)[0]
	return k
	
	
def rsa_enc(p,q, n):

    s = prime_factors(p)
    d = s.find('1*'+str(p))
    if (d==-1):
        return "p is not prime"

    s = prime_factors(q)
    d = s.find('1*'+str(q))
    if (d==-1):
        return "q is not prime"

    if p%3 != 2 or q%3 != 2:
        return "Choose better primes"

    k = p*q
    print "Public key is ", k

    x = n**3 %k
    return "The encrypted number is ", x

def rsa_dec(p,q,t):

    if p%3 != 2 or q%3 != 2:
        return "Choose better primes"

    k = p*q
    s = (2*(p-1)*(q -1)+1)/3
    print s
    c = t**s%k
    return "The decrypted number is ", c


def rsa_pub(public_key, n):

    k = public_key
    print "Public key is ", k

    x = n**3 %k
    return "The encrypted number is ", x


def fact(x): return (1 if x==0 else x*fact(x-1))

def mysin(x,p):
        # compute the sin function using p terms
        s = 0
        j=0
        for i in range(1,p,2):
                k = ((-1.)**j)*x**i /fact(i)
                #print k
                j = j+1
                s = s+k
        return s
                
        
def hw(n = 1):		#It plots n terms of its Fourier Series
	import pylab
	from numpy import pi, sin
	import numpy
	t = numpy.arange(0,2,0.001)
	s = 0
	for i in range(1,n,2):
		s = s + (8/(pi**2))*((-1)**((i-1)/2))*(sin(pi*i*t))/(i**2)
	pylab.grid(True)
	pylab.plot(t,s,lw=3)
	pylab.show()
	



def hw6(n):		#It plots n terms of its Fourier Series
	import pylab
	from numpy import pi, sin
	import numpy
	s = pi**2
	t = 0
	while (t <=2):
		g = 8*sin(pi*t)/s - 8*sin(3*pi*t)/(9*s) + 8*sin(5*pi*t)/(25*s) - 8*sin(7*pi*t)/(49*s)
		print t, g
		t = t + 2./n

def WhatDoIdo(m,n):
        z = 0
        while (m !=0):
                if m%2 == 1:
                        z = z+n
                m = m/2
                n = 2*n
        return z

def fac(n):
        '''
        returns n facorial
        '''
        if n ==0:
                return 1
        else:
                return n*fac(n-1)

        
def sine(x, n):
        '''
        computes the sine of x, where x is in degrees, up to n terms in taylor's seriex
        '''
        from math import pi
        x = 2*pi*x/360
        s = 0
        j = 0
        for i in range(1,n,2):
                h = ((-1)**j)*(x**(i))/fac(i)
                #print i, j, h
                s = s + h
                #print s
                j = j + 1
        return s


def hw8(n = 16):		#2009, hw8
	import pylab
	from numpy import pi, sin
	import numpy
	t = numpy.arange(-2*pi,2*pi,0.001)
	s = 0
	for i in range(1,n):
		s = s + ((-1)**(i-1))*sin(i*t)/i
	pylab.grid(True)
	pylab.plot(t,2*s,lw=3)
	pylab.show()


def hw8v(n,t):		#2009, hw8
	import pylab
	from numpy import pi, sin
	import numpy
	#t = numpy.arange(0,2*pi,0.001)
	s = 0
	for i in range(1,n):
		s = s + ((-1)**(i-1))*sin(i*t)/i
	print 2*s

def lec():		#2009, hw8
	import pylab
	from numpy import pi, sin
	import numpy
	t = numpy.arange(0,1,0.001)
	pylab.grid(True)
	pylab.plot(t,sin(50*t),lw=3)
	pylab.show()

def luhn(n='5444000000000003'):
	'''
	Checks if credit card is valid.
	'''
	n = list(n)
	xf = '2121212121212121'
	xfl = list(xf)
	s = 0
	for i in range(len(n)):
		a = int(n[i])*int(xfl[i])
		if a < 10:
		       s = s + a
		else:
		       s = s + int(str(a)[0]) + int(str(a)[1])
	print 'Total is = ', s
	if s%10 ==0:
		print 'It is a valid card.'
	else:
		print 'The card is not valid.'

def entropy(p=[]):
        '''
        It computes the entropy. 
        '''
        #English p= [0.0800,0.0130,0.0220,0.0460,0.1240,0.0220,0.0200,0.0650,0.0670,0.0010,0.0070,0.0360,0.0250,0.0700,0.0760,0.0160,0.0010,0.0610,0.0620,0.0890,0.0270,0.0080,0.0230,0.0100,0.0200,0.0090]
        #French p = [0.0800,0.0070,0.0350,0.0390,0.1670,0.0120,0.0110,0.0050,0.0760,0.0030,0.0010,0.0490,0.0290,0.0790,0.0580,0.0300,0.0110,0.0740,0.0820,0.0730,0.0550,0.0140,0.0000,0.0060,0.0020,0.0020]
        #German p = [0.0600,0.0170,0.0270,0.0540,0.1800,0.0160,0.0320,0.0410,0.0810,0.0030,0.0130,0.0330,0.0230,0.1060,0.0270,0.0080,0.0000,0.0720,0.0690,0.0570,0.0460,0.0090,0.0150,0.0000,0.0000,0.0110]
        #rolling dice: k = [1,2,3,4,5,6,5,4,3,2,1] and then divide by 36
        #p = [1180./2635, 504./2635, 268./2635, 208./2635, 127./2635, 101./2635, 73./2635, 51./2635, 44./2635, 30./2635, 11./2635, 12./2635, 6./2635, 7./2635, 8./2635, 1./2635, 1./2635, 2./2635, 1./2635]
        from math import log, ceil
        e = 0
        for i in range(len(p)):
                if p[i] !=0:
                        e = e - p[i]*log(p[i],2)
        s = ceil(log(len(p),2))
        print "The entropy is ", e
        print "Symbol length in bits is ", s
        print "The normalized entropy is ", e/s


def bstring(n):
        '''
        It returns a list of all n-bit long binary strings
        '''
        k = ['0','1']
        i = 1
        l = []
        while i < n:
                for j in k:
                        l.append('0'+j)
                for j in k:
                        l.append('1'+j)
                import copy
                k = copy.deepcopy(l)
                l=[]
                i = i+ 1
        return k


def lz(p="iamgladthisisover"):
        '''
        It returns a dictionary words for LZ compression.
        '''
        from math import ceil, log
        #p = "thegovernmentofthepeoplebythepeopleforthepeople"
        #p = "abcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabc"
        #p = "1634733645809253848443133883865090859841783670033092312181110852389333100104508151212118167511579"
        #p = "16347336458092538484431338838650908598417836700330923121811108523893331001045081512121181675115792222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
        #p = "The Democratic effort to secure the two hundred and sixteen votes needed for passage of the legislation came together only after last-minute negotiations involving the White House, the House leadership and a group of Democratic opponents of abortion rights, led by Representative Bart Stupak of Michigan. On Sunday afternoon, members of the group announced that they would support the legislation after Mr  Obama promised to issue an executive order to ensure that federal funds are not used for abortion services."
        #p = "thesoftwaremakerhasstartedgraftingpopularscientificdatabasesandtoolsontoitswindowsazurecloudcomputingservicethisbasicallymeansthatresearchersinvariousfieldsgetaccesstoafastsupercomputeroftheirveryownandcanposequeriestoenormousdatasetsthatmicrosoftkeepsuptodateforthetimebeingmicrosoftwillallowresearchgroupstoperformtheirworkfreeratherthanrentingtimeonazureviaacreditcard"
        L = []            #Comma separated fields

        k = []                  #Distinct symbols
        D=[""]  #Dictionary
        for symbol in p:
                if symbol not in k:
                        k.append(symbol)
        #print k
        ks = bstring(int(ceil(log(len(k),2))))
        #print ks
        kd = {}
        for x in range(len(k)):
                kd[k[x]] = ks[x]
        #print kd
        print '#######################################'
        print "Here is the input:"
        print
        print p
        print
        print "Input length in chars is ", len(p)
        print "There are ", len(k), " distinct symbols in the input."
        t = int(max(1, ceil(log(len(k),2))))*len(p)     #max is used to correct the situation of having one symbol
        print "Input length in bits is", t
        i = 1
        while i <= len(p):
                if p[0:i] in L:
                        i = i + 1
                        continue
                else:
                        L.append(p[0:i])
                        if p[0:i-1] not in  D:
                                D.append(p[0:i-1])
                        p = p[i:]
                        #print p
                        i = 1
        if p != '':
                L.append(p)
        print
        print "Here is the comma separated fields:"
        print L
        print  "comma separated fields length in words is ", len(L)
        print
        print "Here is the Dictionary:"
        print D
        print  "Dictionary length in words is ", len(D)
        print
        Ds = bstring(int(ceil(log(len(D),2))))
        Dd = {}
        for x in range(len(D)):
                Dd[D[x]] = Ds[x]
        print
        c = (int(ceil(log(len(D),2))) + int(ceil(log(len(k),2))))*len(L)
        print "Dictionary encoding requires at most ", c, " bits."
        print "Potential compression saving (not including the header part) in bits is ", t - c
        print "Compression percentage is ", 100*(t-c)/float(t)
        print
        #Final encoding
        E=''
        for word in L:
                #print word
                ws = Dd[word[0:-1]]+ kd[word[-1]]
                #print ws
                E = E + ws + ' '
        print '#######################################'
        print 'symbols table'
        for (x,y) in kd.items():
                print x, "\t", y
        print '#######################################'
        print 'Dictionary table'
        for (x,y) in Dd.items():
                if x == "":
                      print 'null', "\t", "\t", "\t", y
                else:
                        print x, "\t", "\t", "\t", y
        print '#######################################'
        print 'Data compressed'
        print E
		
def sd(n):
	'''
	Rolling two dice n times and computing the probabilty of the sum
	'''
	import random
	p = [1,2,3,4,5,6]
	k = [0,0,0,0,0,0,0,0,0,0, 0, 0, 0]
	for i in range(n):
		f = random.sample(p,1)
		s = random.sample(p,1)
		r = f[0]+s[0]
		k[r] = k[r] + 1
	for i in range(2,13):
		print "Probability of ", repr(i).rjust(2), " is ", k[i]/float(n)
		

def permute(inputData, outputSoFar):
	for elem in inputData:
		if elem not in outputSoFar:
			outputSoFar.append(elem)
			if len(outputSoFar) == len(inputData):
				print outputSoFar
			else:
				permute(inputData, outputSoFar) # --- Recursion
			outputSoFar.pop() 

          
                        
                        
class Permutation:
	def __init__(self, justalist):
		self._data = justalist[:]
		self._sofar = []
	def __iter__(self):
			return self.next()
	def next(self):
		for elem in self._data:
			if elem not in self._sofar:
				self._sofar.append(elem)
				if len(self._sofar) == len(self._data):
					yield self._sofar[:]
				else:
					for v in self.next():
						yield v
				self._sofar.pop()  
        
        
def h8y11():
        a = '141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609'
        k = list(a)
        print len(k)
        l= []
        for i in range(10):
                l.append(k.count(str(i)))
        print l
        p = []
        for i in l:
                p.append(round(float(i)/len(k), 3))
        print p
        h = 0.
        from math import log
        for i in p:
                h = h + (-i*log(i,2))
        print 'Entropy is ', h
        print 'Normalized Entropy is ', h/4.
        lz(a)
        
        
        
def Karprekar():
        '''
        The mysterious 6174
        '''
        s = raw_input('Enter a four digit number not like 1111, 2222, ... = ')
        i = 0
        if sorted(s) == sorted(s, reverse = True):
                print 'Bad input'
                return
        while True:
                h = ''.join(sorted(s, reverse = True))
                l = ''.join(sorted(s))
                d = int(h) - int(l)
                print int(h), ' - ', int(l), ' = ', d
                if str(d) == s:
                        print s, i
                        break
                else:
                        s = str(d)
                        i = i+1
                        #raw_input()

def randPwd(length = 5):
        import string
        import random
        
        UC = string.ascii_uppercase
        pwd=''
        for i in range(length):
                pwd = pwd + random.choice(UC)
        return pwd

def randPixel():
        import random
        r = random.choice(range(256))
        g = random.choice(range(256))
        b = random.choice(range(256))
        return((r,g,b))

def rgbFlag():
        import numpy as np
        picture = np.empty( (300,400), dtype=object)
        for i in range(100):
                for j in range(400):
                        picture[i,j] = (255,0,0)
        for i in range(100, 200):
                for j in range(400):
                        picture[i,j] = (0,255,0)
        for i in range(200, 300):
                for j in range(400):
                        picture[i,j] = (0,0,255)                       
        return picture

def randPicture(row = 3, column = 4):
        import numpy as np
        picture = np.empty( (row,column), dtype=object)
        for i in range(row):
                for j in range(column):
                        picture[i,j] = randPixel()
        return picture

def picSave(picture=rgbFlag()):
        '''
        Saves picture in BMP format
        '''
        import Image, ImageDraw
        import numpy as np
        
        row,column = picture.shape
        #PIL picture format; one "plane" for each color
        pixels = np.zeros((row,column,3), dtype=np.uint8)
        #print pixels
        for i in range(row):
                for j in range(column):
                        pixels[i,j,0] = picture[i,j][0]
                        pixels[i,j,1] = picture[i,j][1]
                        pixels[i,j,2] = picture[i,j][2]
        
        im = Image.fromarray(pixels)
        filePath = raw_input('Enter file path: ')
        out = open(filePath + '\myPic.bmp', 'w')
        im.save(out, "BMP")
        out.close()
        print 'Your picture is stored in: ', filePath + '\myPic.bmp'
  

	
