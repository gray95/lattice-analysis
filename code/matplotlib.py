# matplotlib plots

def plot_corr( corr):   
    print ("Reading data from  " , filename_in)
    data = open(filename_in, 'r')
    
    corr = np.genfromtxt(data, autostrip=True)
    
    corr_mean = np.mean(corr, axis=0)
    uncertainty = np.std(corr, axis=0)
    
    corr_mean = corr_mean[1:]
    uncertainty = uncertainty[1:]
    
    no_config = corr.shape[0]
    nt        = corr.shape[1]
    
    tt, meff, meff_err = calc_meff(corr, nt, no_config) 
    
    ### plot the data
    
    plt.subplot(311)
    plt.errorbar(range(len(corr_mean)),[abs(corr_mean[q]) for q in range(len(corr_mean))], fmt='bo', yerr=None  , alpha=0.5, markersize=2)
    plt.yscale('log', nonposy='clip')
    plt.title(r'Mean 1-+ Correlator')  
    plt.ylabel(r'$\langle 0|e^{-HT}|0 \rangle$')
    
    E = []
    for i in range(len(corr_mean)-1):
        hldr = m.log(abs((corr_mean[i]/corr_mean[i+1])))
        E.append(hldr)
        
    plt.subplot(312)
    plt.plot(range(len(E)), [E[x] for x in range(len(E))], 'b+')
    plt.ylabel(r'$\log( \frac{G(t)}{G(t+a)})$')
    #plt.xlabel('time')
    plt.axis([0,47,-6,6])
    
    plt.subplot(313)
    plt.errorbar(tt, meff , meff_err  ,  fmt= 'ro', markersize=2, alpha=0.5)
    plt.xlabel('t')
    plt.ylabel('meff')
    plt.xlim(0,9)
    plt.ylim(-2,6)
    
    if SAVEFIG :
        plt.savefig('test')
    
    
    plt.show()
        
        