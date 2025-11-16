#Reimportar las librerías necesarias
import numpy as np
import pandas as pd
import math as mt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
import os
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
#DATOS EXPERIMENTALES
def extraer_vectores(nombre_archivo, delimitador='\t'):
    """
    Lee un archivo de texto y devuelve los vectores H y B
    ignorando columnas adicionales: coge las 2 primeras.
    """
    # Lee el archivo omitiendo encabezados no numéricos
    df = pd.read_csv(nombre_archivo, sep=delimitador, skiprows=2, engine='python')
    #Toma solo las dos primeras columnas
    df = df.iloc[:, :2]
    df.columns = ['H', 'B']
    #Conversion punto y float
    df['H'] = df['H'].astype(str).str.replace(',', '.').astype(float)
    df['B'] = df['B'].astype(str).str.replace(',', '.').astype(float)
    return df

#CÁLCULOS
"""
1-Gráficar M(H)=>uM=B/uH
H=I*2586
¿Se alcanza saturación?
"""
def ciclosimple():
    def ciclo(df0, nombre='M(H)',nciclo='Ciclo N'):
        H=df0.iloc[:, 0].to_numpy()
        H=(2586/9524)*H
        B=df0.iloc[:, 1].to_numpy()
        mu0=4*mt.pi*1e-7
        M=(B/mu0)-H
        I=H/9524
        X=H
        Y=M
        plt.scatter(X, Y, label="Datos experimentales", color="red", s=1.5)
        plt.xlabel("H (A/m) ")
        plt.ylabel("M (A/M)")
        plt.title(f"M en función de H ({nciclo})")
        plt.legend()
        archivo_png = f"{nombre}.png"
        plt.savefig(archivo_png, dpi=300)
        print(f"✅ Gráfico guardado como: {archivo_png}")
        data={
            'I (A)':I,
            'B (T)':B,
            'M (A/M)':M,
            'H (A/M)':H
        }
        plt.close()
        df=pd.DataFrame(data)
        return df
    df1=extraer_vectores(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo1bueno.txt")
    dfC1=ciclo(df1,'M1(H1)','ciclo 1')
    df2=extraer_vectores(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo2bueno.txt")
    dfC2=ciclo(df2,'M2(H2)','ciclo 2')
    df5=extraer_vectores(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclos de 5.txt")
    return dfC1, dfC2, df5
def editar_puntos_txt(ruta_txt, tol=0.02, save_path=None):
    """
    Permite visualizar y editar (eliminar) pares de datos de las dos primeras columnas de un archivo .txt.
    - ruta_txt: ruta al archivo .txt (debe tener al menos dos columnas numéricas).
    - tol: tolerancia para asociar un click a un punto (0..1, relativo).
    - save_path: ruta donde guardar el archivo modificado (opcional, si None pregunta).
    """
    # Leer solo las dos primeras columnas numéricas
    df= pd.read_csv(ruta_txt, sep='\t', skiprows=2, engine='python', header=None)
    df= df.iloc[:, :2]
    df.columns=['Col1', 'Col2']
    X=df['Col1'].astype(str).str.replace(',', '.').astype(float)
    Y=df['Col2'].astype(str).str.replace(',', '.').astype(float)
    range_X = np.ptp(X) if np.ptp(X) != 0 else 1.0
    range_Y = np.ptp(Y) if np.ptp(Y) != 0 else 1.0
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(X, Y, c='red', s=8)
    ax.plot(X,Y,c='green', linewidth=1.2)
    ax.set_xlabel("Columna 1")
    ax.set_ylabel("Columna 2")
    ax.set_title("Click sobre puntos a eliminar. Presiona Enter para terminar")
    plt.draw()
    try:
        clicks=plt.ginput(n=0, timeout=0)
    except Exception as e:
        print(f"Error interactivo: {e}")
        plt.close(fig)
        return
    if not clicks:
        print("Ningún punto seleccionado. No se realizaron cambios.")
        plt.close(fig)
        return
    X_min,Y_min=X.min(),Y.min()
    X_norm=(X-X_min)/range_X
    Y_norm=(Y-Y_min)/range_Y
    removed_idx = set()
    for (x_click, y_click) in clicks:
        x_n=(x_click-X_min)/range_X
        y_n=(y_click-Y_min)/range_Y
        dists=np.sqrt((X_norm-x_n)**2+(Y_norm-y_n)**2)
        idx=int(np.argmin(dists))
        if dists[idx] <= tol:
            removed_idx.add(idx)
            ax.plot(X.iloc[idx], Y.iloc[idx], marker='x', color='black', markersize=8)
            plt.draw()
    if not removed_idx:
        print("Ningún click se asoció a puntos (ajusta tolerancia si es necesario).")
        plt.show()
        plt.close(fig)
        return
    df_clean=df.drop(index=list(removed_idx)).reset_index(drop=True)
    print(f"Se eliminaron {len(removed_idx)} puntos: índices {sorted(list(removed_idx))}")
    if save_path is None:
        respuesta=input(f"¿Desea sobrescribir el archivo original '{ruta_txt}'? (s/n): ").strip().lower()
        if respuesta == "s":
            save_path=ruta_txt
        else:
            save_path=None

    if save_path:
        try:
            df_clean.to_csv(save_path, sep='\t', index=False, header=False, float_format='%.6f')
            print(f"✅ Datos guardados en: {save_path}")
        except Exception as e:
            print(f"⚠️ Error al guardar: {e}")
    plt.show()
    plt.close(fig)
    return df_clean
"""
2-Parámetros magneticos
    a)Campo coercitivo
    b)Saturación magnética y remanencia
    c)Permeabilidad
"""
def analisisciclo(dfC1, dfC2):
    def HcMsMr(df,ciclo="ciclo 0"):
        df = df.copy()
        df["M"] = df.iloc[:, 2]
        df["H"] = df.iloc[:, 3]
        M_all = df["M"].to_numpy()
        H_all = df["H"].to_numpy()
        dh=np.diff(H_all)
        dh = np.convolve(dh, np.ones(3)/3, mode='same')
        sign_changes=np.where(np.sign(dh[:-1]) * np.sign(dh[1:]) < 0)[0]
        if sign_changes.size > 0:
            primer_max_idx = sign_changes[0] + 1
        else:
            primer_max_idx = int(np.argmax(np.abs(H_all)))
        df = df.loc[primer_max_idx:].reset_index(drop=True)
        M = df["M"].to_numpy()
        H = df["H"].to_numpy()
        def HC(M,H):
            def interpol(x,y,target=0):
                xs=[]
                for i in range(len(x)-1):
                    y1,y2=y[i],y[i+1]
                    if y1==target:
                       xs.append(x[i])
                    if (y1*y2<0) or (y2==target):
                        t = (target - y1) / (y2 - y1)
                        x_interp = x[i] + t*(x[i+1]-x[i])
                        xs.append(x_interp)
                return np.array(xs, dtype=float)
            Hc=interpol(H,M,0)
            Hcmean=np.mean(np.abs(Hc))
            #print(7205*0.521965)
            Hcmean=Hcmean*4*np.pi*1e-7
            if Hcmean<1e-3:
                Material=f"BLANDO"
            if Hcmean>0.1:
                Material=f"DURO"
            return Hc, Hcmean, Material
        def MsMr(df):
            def saturado(df):
                df_sortM=df.sort_values(by="M")
                Saturadoneg=df_sortM.head(10)
                Msn=Saturadoneg["M"].mean()
                Hsn=Saturadoneg["H"].mean()
                Saturadopos=df_sortM.tail(10)
                Msp=Saturadopos["M"].mean()
                Hsp=Saturadopos["H"].mean()
                Satnp=np.array([Msp,Msn,Hsp,Hsn])
                Sat=np.array([((Msp+abs(Msn))/2,(Hsp+abs(Hsn))/2)])
                return Satnp, Sat
            def remanente(df):
                # interpola M en H=0 para todos los cruces y devuelve un único Mr positivo y uno negativo
                H = df.iloc[:, 3].to_numpy(dtype=float)
                M = df.iloc[:, 2].to_numpy(dtype=float)
                mr_list = []
                for i in range(len(H) - 1):
                    H1, H2 = H[i], H[i+1]
                    M1, M2 = M[i], M[i+1]
                    # caso exacto en un punto
                    if H1 == 0:
                        mr_list.append(M1)
                    # cruce entre H1 y H2
                    if H1 * H2 < 0:
                        # evitar división por cero
                        if (H2 - H1) != 0:
                            t = -H1 / (H2 - H1)
                            M_interp = M1 + t * (M2 - M1)
                            mr_list.append(M_interp)
                    # considerar H2==0 para no perder último punto
                    if H2 == 0:
                        mr_list.append(M2)
                if len(mr_list) == 0:
                    pos_val = 0.0
                    neg_val = 0.0
                else:
                    mr_arr = np.array(mr_list, dtype=float)
                    pos_candidates = mr_arr[mr_arr > 0]
                    neg_candidates = mr_arr[mr_arr < 0]
                    pos_val = float(pos_candidates.max()) if pos_candidates.size > 0 else 0.0
                    neg_val = float(neg_candidates.min()) if neg_candidates.size > 0 else 0.0
                Mr = np.array([pos_val, neg_val], dtype=float)  # [Mr_pos, Mr_neg]
                Mrmean = float(np.mean(np.abs(Mr)))
                return Mr, Mrmean
            Satnp, Sat=saturado(df)
            Mr, Mrmean=remanente(df)
            Mr = np.array(Mr, dtype=float)
            return Satnp, Sat, Mr, Mrmean
        Hc, Hcmean, material=HC(M,H)
        Satnp, Sat, Mr, Mrmean=MsMr(df)
        
        plt.figure(figsize=(10, 8))
        plt.scatter(H, M, label="Datos experimentales", color="red", s=1.5)
        plt.title(f"Parámetros magnéticos del {ciclo}")
        # Dibujar ejes en x=0 e y=0
        plt.axhline(y=0, linestyle='--', color='black', linewidth=0.8)
        plt.axvline(x=0, linestyle='--', color='black', linewidth=0.8)

        # Marcar +/-Hc
        for hc_val in Hc:
            plt.plot(hc_val, 0, 'o', color='aquamarine', markersize=8)
        plt.plot([], [], 'o', color='aquamarine', markersize=8, label="Hc")
        
        # Marcar +/-Mr
        for mr_val in Mr:
            plt.plot(0, mr_val, 's', color='green', markersize=8)
        plt.plot([], [], 's', color='green', markersize=8, label="M remanente")
        
        # Marcar +/-Ms
        plt.plot(Satnp[2], Satnp[0], 'o', color='purple', markersize=8)
        plt.plot(Satnp[3], Satnp[1], 'o', color='purple', markersize=8)
        plt.plot([], [], 'o', color='purple', markersize=8, label="M saturación")
        mu0=4*np.pi*1e-7
        parametrosmag={
            "Hc+/- [A/m]":Hc,
            "Hcmean[A/m]":Hcmean/mu0,
            "Hcmean[T]":Hcmean,
            "Tipo material":material,
            "(+/-mu0Ms,+/-mu0Hs) [T]":Satnp*mu0,
            "(mu0Ms,mu0Hs) [T]":Sat*mu0,
            "(+/-Ms,+/-Hs) [A/m]":Satnp,
            "(Ms,Mr) [A/m]":Sat,
            "+/-mu0Mr [T]":Mr*mu0,
            "mu0Mr [T]":Mrmean*mu0,
            "+/- Mr [A/m]":Mr,
            "Mr [A/m]":Mrmean
        }
        
        plt.xlabel("H (A/m)")
        plt.ylabel("M (A/m)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(f"Parametros_magneticos_{ciclo}.png", dpi=300)
        plt.close()
        print(f"✅ Gráfico guardado como: Parametros magneticos_{ciclo}.png")
        dfpm=pd.DataFrame([parametrosmag])        
        return dfpm
    def CurvIni(df,puntos=232, derivada=50,ciclo="CICLO_0"):
        #Seleccion a mano de la curva de imanación inicial (232 puntos C1, 234 en C2)
        df_ini=df.iloc[:puntos].copy()
        H=df_ini.iloc[:, 3].to_numpy(dtype=float) 
        M=df_ini.iloc[:, 2].to_numpy(dtype=float) 
        # G1:Curva inicial de imanación
        plt.figure(figsize=(8, 6))
        plt.scatter(H, M, color='b', label=f'Datos experimentales ({puntos} pts)', s=8)
        #Ajuste pol grado: 5
        grado=5
        coef=np.polyfit(H, M, grado)
        p=np.poly1d(coef)
        #Curva de ajuste suave
        H_fit=np.linspace(H.min(), H.max(), 300)
        M_fit=p(H_fit)
        plt.plot(H_fit, M_fit, 'r-', label=f'Ajuste polinómico (grado {grado})')
        plt.xlabel("H (A/m)")
        plt.ylabel("M (A/m)")
        plt.title(f"Curva inicial de imanación ({ciclo})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"Curva_Inicial_{puntos}pts_{ciclo}.png", dpi=300)
        plt.close()
        print(f"✅ Gráfico guardado como: Curva_Inicial_{puntos}pts_{ciclo}.png")
        #Se eligen los puntos de variable=derivada equiespaciados para derivar
        H_50=np.linspace(H.min(), H.max(), derivada)
        M_50=p(H_50)
        #Se deriva
        p_deriv=np.polyder(p)
        dM_dH_50=p_deriv(H_50)

        #G2: Susceptibilidad magnética diferencial
        plt.figure(figsize=(8, 6))
        plt.plot(H_50, dM_dH_50, 'go', markersize=3, label='Susceptibilidad magnética (dM/dH)')
        plt.xlabel("H (A/m)")
        plt.ylabel("χ=dM/dH")
        plt.title(f"Susceptibilidad magnética respecto al Campo H ({ciclo})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"Susceptibilidad_vs_H_{derivada}pts_{ciclo}.png", dpi=300)
        plt.close()
        print(f"✅ Gráfico guardado como: Susceptibilidad_vs_H_{derivada}pts_{ciclo}.png")

        #Guardar datos en dfresult
        dfresult = pd.DataFrame({
            "H (A/m)":H_50,
            "M_ajuste (A/m)":M_50,
            "dM/dH":dM_dH_50
        })
        return dfresult
    dfpm1 = HcMsMr(dfC1,"ciclo 1")
    dfpm2 = HcMsMr(dfC2,"ciclo 2")
    dfresult1 = CurvIni(dfC1,232,50,"ciclo 1")
    dfresult2 = CurvIni(dfC2,234,50,"ciclo 2")
    return dfpm1, dfpm2, dfresult1, dfresult2
"""
4-Hacer 5 ciclos y representar los máximos
"""
def cinco_ciclos(df):
    #Procesar el DataFrame df0 para obtener los resultados de los 5 ciclos
    def procesinutil(df):
        df0 = df.copy()
        H=df0.iloc[:, 0].to_numpy()
        H=(2586/952.4)*H
        B=df0.iloc[:, 1].to_numpy()
        mu0=4*mt.pi*1e-7
        M=(B/mu0)-H
        I=H/952.4
        return pd.DataFrame({'H (A/m)': H, 'M (A/m)': M, 'I (A)': I, 'B (T)': B})
    def selec_ciclos(df,n0,n1):
        df0 = df.iloc[n0:n1].copy()
        H=df0.iloc[:, 0].to_numpy()
        H=(2586/952.4)*H
        B=df0.iloc[:, 1].to_numpy()
        mu0=4*mt.pi*1e-7
        M=(B/mu0)-H
        I=H/952.4
        dfr=pd.DataFrame({'H (A/m)': H, 'M (A/m)': M, 'I (A)': I, 'B (T)': B})
        return dfr
    def maximos(dfr):
        dfr_sortM=dfr.sort_values(by="M (A/m)")
        Saturado=dfr_sortM.tail(5)
        Ms=Saturado['M (A/m)'].mean()
        Hs=Saturado['H (A/m)'].mean()
        return np.array([Hs,Ms])
    def CII(dfr, dfsat):
        H=dfr['H (A/m)'].to_numpy()
        M=dfr['M (A/m)'].to_numpy()
        Hs=dfsat['Hs (A/m)'].to_numpy()
        Ms=dfsat['Ms (A/m)'].to_numpy()
        plt.scatter(H, M, color='blue', label='Curva de imanación inicial', s=1.5)
        plt.scatter(Hs, Ms, color='black', label='Puntos máximos')
        plt.xlabel("H (A/m)")
        plt.ylabel("M (A/m)")
        plt.legend()
        plt.title("Curva de imanación inicial")
        plt.grid(True, alpha=0.3)
        plt.savefig(f"Curva_Imanacion_Inicial_MAX.png", dpi=300)
        plt.show()
        print(f"✅ Gráfico guardado como: Curva_Imanacion_Inicial_MAX.png")
        return
        
    dfII =selec_ciclos(df,0,357)
    dfr1 = selec_ciclos(df,357,1163)
    Sat1 = maximos(dfr1)
    dfr2 = selec_ciclos(df,1163,1970)
    Sat2 = maximos(dfr2)
    dfr3 = selec_ciclos(df,1970,2573)
    Sat3 = maximos(dfr3)
    dfr4 = selec_ciclos(df,2573,3027)
    Sat4 = maximos(dfr4)
    dfr5 = selec_ciclos(df,3027,len(df))
    Sat5 = maximos(dfr5)
     #Procesar los 5 ciclos 
    ciclos = [
        (357, 1163),
        (1163, 1970),
        (1970, 2573),
        (2573, 3027),
        (3027, len(df))
    ]

    colores = ['r', 'g', 'b', 'orange', 'purple']
    sats = []

    plt.figure(figsize=(10, 7))

    for i, (n0, n1) in enumerate(ciclos, start=1):
        dfr = selec_ciclos(df, n0, n1)
        Sat = maximos(dfr)
        sats.append(Sat)

        # Graficar el ciclo
        plt.scatter(dfr['H (A/m)'], dfr['M (A/m)'],
                 color=colores[i - 1], s=1.5,
                 label=f'Ciclo {i}')

        # Punto de saturación
        plt.plot(Sat[0], Sat[1], 'ko', markersize=6)
        offset_x = (dfr['H (A/m)'].max() - dfr['H (A/m)'].min()) * 0.01
        offset_y = (dfr['M (A/m)'].max() - dfr['M (A/m)'].min()) * 0.015
        plt.text(Sat[0] + offset_x, Sat[1] + offset_y,
                 f"Ms={Sat[1]:.2e}", fontsize=8, color='black', ha='left', va='bottom')
    plt.xlabel("H (A/m)")
    plt.ylabel("M (A/m)")
    plt.title("Curvas de magnetización - 5 ciclos")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Ciclos_5_Saturacion.png", dpi=300)
    plt.show()

    print("✅ Gráfico guardado como: Ciclos_5_Saturacion.png")

    #Retornar datos max
    sats_df = pd.DataFrame(sats, columns=["Hs (A/m)", "Ms (A/m)"])
    sats_df.index = [f'Ciclo {i}' for i in range(1, 6)]
    #print(sats_df)
    CII(dfII, sats_df)
    inutil=procesinutil(df)
    return sats_df,inutil

##CALL
#editar_puntos_txt(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo1bueno.txt", tol=0.02, save_path=None)
#editar_puntos_txt(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo2bueno.txt", tol=0.02, save_path=None)
#editar_puntos_txt(r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclos de 5.txt",tol=0.02, save_path=None)
dfC1, dfC2, df5 = ciclosimple()
dfpm1, dfpm2, dfresult1, dfresult2 =analisisciclo(dfC1, dfC2)
sats_df, inutil = cinco_ciclos(df5)
#print('\nCICLO 1\n',dfpm1,'\nCICLO 2\n',dfpm2)
# Lista de archivos a procesar
archivos = [
    r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclos de 5.txt",
    r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo 5V.txt",
    r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo2bueno.txt",
    r"C:\Users\JuanR\Desktop\FÍSICA\3º\LAB FIS III\ELECTRO\HISTERESIS\ciclo1bueno.txt"
]

# Diccionario para el JSON
datos_json = {}

# Crear Excel con varias hojas
with pd.ExcelWriter('vectores_H_B.xlsx') as writer:
    for archivo in archivos:
        try:
            nombre_simple = archivo.split('\\')[-1].replace('.txt', '')
            df = extraer_vectores(archivo)
            df.to_excel(writer, sheet_name=nombre_simple, index=False)
            datos_json[nombre_simple] = {'H': df['H'].tolist(), 'B': df['B'].tolist()}
            print(f"✅ {nombre_simple} procesado correctamente.")
        except FileNotFoundError:
            print(f"⚠️ Archivo no encontrado: {archivo}")
        try:    
                dfC1.to_excel(writer, sheet_name="Ciclo1", index=False)
                dfC2.to_excel(writer, sheet_name="Ciclo2", index=False)
                inutil.to_excel(writer, sheet_name="Ciclos5", index=False)
                dfpm1.to_excel(writer, sheet_name="Parametros_Ciclo1", index=False)
                dfpm2.to_excel(writer, sheet_name="Parametros_Ciclo2", index=False)
                dfresult1.to_excel(writer, sheet_name="ResultadosX_Ciclo1", index=False)
                dfresult2.to_excel(writer, sheet_name="ResultadosX_Ciclo2", index=False)
                sats_df.to_excel(writer, sheet_name="Maximos_CincoCiclos", index=False)
        except Exception as e:
            print(f"⚠️ Error al procesar analisisciclo: {e}")
# Guardar todo en un JSON
with open('vectores_H_B.json', 'w', encoding='utf-8') as f:
    json.dump(datos_json, f, ensure_ascii=False, indent=4)

# Guardar ambos DataFrames en nuevas hojas del mismo Excel

print("✅ Vectores H y B guardados en:")
print("  - Excel: vectores_H_B.xlsx")
print("  - JSON:  vectores_H_B.json")