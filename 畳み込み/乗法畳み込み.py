"""乗法畳み込み C_k = Σ_{ij≡k (mod p)} A_i B_j。

原始根で指数に直し、NTT 畳み込みする。

使用例:
    mul_convolution([1, 2], [3, 4], 5, 998244353)
"""
def convolution(f,g,mod):
    
    """NTT + Garner による任意法畳み込み。"""
    def _convolution(f,g,_mod):

        """1 つの NTT 素数 _mod で畳み込む。"""
        n = len(bin(len(f)+len(g)-1)) - 2
        fft_length = 1<<n

        f = f + [0]*(fft_length - len(f))
        g = g + [0]*(fft_length - len(g))

        if _mod==998244353:
            w = [1, 998244352, 911660635, 372528824, 929031873, 452798380, 922799308, 781712469, 476477967, 166035806, 258648936, 584193783, 63912897, 350007156, 666702199, 968855178, 629671588, 24514907, 996173970, 363395222, 565042129, 733596141, 267099868, 15311432]
            iw = [1, 998244352, 86583718, 509520358, 337190230, 87557064, 609441965, 135236158, 304459705, 685443576, 381598368, 335559352, 129292727, 358024708, 814576206, 708402881, 283043518, 3707709, 121392023, 704923114, 950391366, 428961804, 382752275, 469870224]
        
        if _mod==897581057:
            w = [1, 897581056, 200527991, 850960045, 227655573, 685177417, 661961559, 717889083, 688546301, 64431346, 762907769, 781659575, 604882016, 471181658, 242773703, 313099125, 288794207, 732004569, 437566725, 430897771, 279727937, 91704119, 523358721, 872686320]
            iw =[1, 897581056, 697053066, 442470459, 502631723, 192192108, 366473218, 285218810, 627498913, 632928577, 715124372, 829482092, 895669752, 835819291, 210274124, 7242324, 530138839, 365592405, 712687518, 812501856, 244025573, 353112847, 793229247, 354917575]

        if _mod==880803841:
            w = [1, 880803840, 121444121, 547680885, 836988352, 170630252, 547743738, 390590270, 755881750, 119481987, 622213777, 634844223, 496183605, 872875137, 41469254, 551868471, 219288049, 198000217, 579409128, 733691905, 566136041, 374515633, 402082372, 273508579]
            iw = [1, 880803840, 759359720, 339414624, 282082127, 83908436, 623501316, 879302015, 26105166, 708522529, 769895303, 843755407, 710708181, 623500536, 528308065, 542164623, 817679620, 571049407, 409417309, 504998132, 352282463, 252040680, 400443141, 109748732]


        def fft(a):

            """数論変換(破壊的)。"""
            for i in range(fft_length):
                j = 0
                for k in range(n):
                    j |= ((i>>k)&1) << (n - 1 - k)

                if i<j:
                    a[i],a[j] = a[j],a[i]
            
            for nn in range(n):
                b = 1<<nn
                wj = 1
                _wj = iw[nn+1]
                for j in range(b):


                    for k in range(0,fft_length,2*b):

                        s = a[j+k]
                        t = a[j+k+b]*wj%_mod
                        a[j+k] = (s+t)%_mod
                        a[j+k+b] = (s-t)%_mod

                    wj = wj*_wj%_mod            
            return a
        
        def ifft(a):
            """逆数論変換。"""
            for i in range(fft_length):
                j = 0
                for k in range(n):
                    j |= ((i>>k)&1) << (n - 1 - k)

                if i<j:
                    a[i],a[j] = a[j],a[i]
            
            for nn in range(n):
                b = 1<<nn
                wj = 1
                _wj = w[nn+1]
                for j in range(b):

                    for k in range(0,fft_length,2*b):

                        s = a[j+k]
                        t = a[j+k+b]*wj%_mod
                        a[j+k] = (s+t)%_mod
                        a[j+k+b] = (s-t)%_mod

                    wj = wj*_wj%_mod            
            inv = pow(fft_length,_mod-2,_mod)
            return [i*inv%_mod for i in a]
                
        F = fft(f)
        G = fft(g)
        H = [i*j%_mod for i,j in zip(F,G)]
        return [i for i in ifft(H)]
    
    f = [i%mod for i in f]
    g = [i%mod for i in g]
    
        
    x = _convolution(f,g,998244353)
    if mod==998244353:
        return x
    y = _convolution(f,g,897581057)
    z = _convolution(f,g,880803841)

    m1 = 998244353
    m2 = 897581057
    m3 = 880803841

    m1_inv_m2 = pow(m1,m2-2,m2)
    m12_inv_m3 = pow(m1*m2,m3-2,m3)
    m12_mod = m1*m2%mod

    res = [0]*len(x)

    for i in range(len(x)):

        v1 = (y[i]-x[i])*m1_inv_m2%m2
        v2 = (z[i]-(x[i]+m1*v1)%m3)*m12_inv_m3%m3
        c3 = (x[i]+ m1*v1 + m12_mod*v2)%mod
        res[i] = c3
    
    return res


def mul_convolution(f,g,p,mod):
    
    #C_k = (sum(A_i*B_j)(k=i*j%p))%modを求める
    # f = [A1,A2,,,,A(p-1)]
    # g = [B1,B2,,,,B(p-1)]
    """C_k = Σ_{i*j ≡ k (mod p)} f[i-1] g[j-1] を mod で返す。"""
    if p==2:
        return [0,f[0]*g[0]%mod]
    #p-1の約数列挙
    _p = p-1
    pp = []
    for j in range(2,int(_p**(1/2)) + 3):
        
        if j**2>p:
            break
        
        if _p%j==0:
            if j**2!=_p:
                pp.append(j)
            pp.append(_p//j)
    #2から順に原始根になるか計算
    for i in range(2,p):
        
        if all(pow(i,s,p)!=1 for s in pp):
            
            p_root = i
            break

    R = [1]*(p)
    for i in range(1,p):
        R[i] = p_root*R[i-1]%p
    R_inv = {r:i for i,r in enumerate(R)}
    F = [0]*p
    G = [0]*p
    for i in range(p-1):
        gg = g[i]
        ff = f[i]
        G[R_inv[i+1]] = gg
        F[R_inv[i+1]] = ff
    h = convolution(F,G,mod)
    H = [0]*p
    for i in range(len(h)):
        hh = h[i]
        idx = ((i-1)%(p-1))+1
        ridx = R[idx]
        H[ridx] = (H[ridx] + hh)%mod

    return H
