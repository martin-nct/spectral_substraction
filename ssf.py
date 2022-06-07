def ssf(x, fs, largo=0.01, solap=None, silencio=0.3, b=1):
    
    '''
    Aplica reducción de ruido por el método de substracción espectral. Requiere que la señal no tenga información
    útil en el inicio. Procesa la señal con reducción de ruido residual.
    
    Entrada:
    x: array 
        Señal a filtrar
    fs: int 
        Frecuencia de muestreo
    largo: float
        Longitud en segundos de la ventana para realizar la STFT. Por defecto largo = 0.01 seg
    solap: float
        El solapamiento en segundos entre ventana y ventana en la STFT. Por defecto es la mitad del largo. Debe ser menor
        al largo.
    silencio: float
        Tiempo en segundos al inicio de la señal donde se estimará el ruido. Se requiere que no haya información útil.
    b: int/float
        Parámetro de sobre-sustracción. Usualmente 1<b<2
    
    Salida:
    salida: array
        Señal filtrada
    tiempo: array
        Vector de tiempo de la señal
    '''
    
    if solap != None:
        solap = int(solap*fs)    # Pone el paso en muestras
    f, t, X = sig.stft(x, fs, nperseg=int(largo*fs), noverlap=solap)    # Calcula STFT
    X = X.T    # Trasponer la matriz
    
    ventanas = []
    fases = []
    for i in range(len(X)):
        ventanas.append(np.abs(X[i]))
        fases.append(np.angle(X[i]))
    
    ruido = np.array(ventanas[0])
    muestras = int(silencio * fs)    # Cantidad de muestras sin señal útil
#     nvent = np.shape(ventanas)[0]    # Número de ventanas totales
#     M = int(muestras/nvent)    # Las primeras M ventanas con ruido
    M = int(muestras/len(ventanas[0]))
    
    for i in range(1,M):
        ruido += ventanas[i]
    ruido /= M    # Estimador de ruido
    
    ruido *= b    # Oversubstract
    
    for i in range(len(ventanas)):
        for k in range(len(ventanas[i])):
            if ventanas[i][k] > ruido[k]:
                ventanas[i][k] -= ruido[k]
            else: ventanas[i][k] = 0
    
    # Reducción de ruido residual:
    
    maxres = ventanas[0]            
    for i in range(1,M):
        for k in range(len(ventanas[i])):
            if ventanas[i][k] > maxres[k]:
                maxres[k] = ventanas[i][k]    # Selecciona las magnitudes máximas en el silencio inicial luego de la substracción
    
    for i in range(1,len(ventanas)-1):
        for k in range(len(ventanas[i])):
            if ventanas[i][k] < maxres[k]:
                ventanas[i][k] = min(ventanas[i+1][k], ventanas[i][k], ventanas[i-1][k])
                # Elige el valor mínimo en ventanas adyacentes para esa frecuencia
    
    
    out = []
    for i in range(len(ventanas)):
        out.append(ventanas[i] * np.exp(1j*fases[i]))    # Agregamos la magnitud y la fase
    out = np.asarray(out)
    
    tiempo, salida = sig.istft(out.T,fs, noverlap=solap)     # Antitransformamos
    return salida, tiempo
