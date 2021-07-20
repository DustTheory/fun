import math
from IPython.core.display import display, HTML
import colorsys

passCount = failCount = totalCount = 0

def resetTests():
    global passCount, failCount, totalCount
    passCount = failCount = totalCount = 0

def displayTotal():
    outputColor = numberToColor(1);
    if(totalCount > 0):   
        outputColor = numberToColor(passCount/totalCount)

    html = f"""
    <div style=\"color:{outputColor};\">
        <b> TESTS SUMMARY </b>
        <br>
        <span style=\"margin-left: 50px;\"> TOTAL: {totalCount} </span> <br>
        <span style=\"margin-left: 50px;\"> PASSED: {passCount} </span> <br>
        <span style=\"margin-left: 50px;\"> FAILED: {failCount} </span> <br>
    </div>
    """

    display(HTML(html))
    

def EvalTest(expected, output, tolerance, description):
    global totalCount, passCount, failCount
    totalCount += 1
    error = 0 if expected == output else abs((expected - output)/expected)
    status = "PASS" if error < tolerance else "FAIL"
    if(status == "PASS"):
        passCount += 1
    else:
        failCount += 1

    outputColor = generateOutputColor(error, tolerance)

    html = f"""
    <div style=\"color:{outputColor};\">
        <b>{status} TEST {totalCount} ({description}): </b>
        <br>
        <span style=\"margin-left: 50px;\">Error: {round(error, 6):2e}, Tolerance: {round(tolerance, 6):2e}</span>
    </div>
    """

    display(HTML(html))


def generateOutputColor(error, tolerance):
    normalizedErr = 1 - max(0, min(1, error / tolerance))
    return numberToColor(normalizedErr)

def numberToColor(i):
    rgb = colorsys.hls_to_rgb(i*0.333, .3, 1) 
    return f"rgb({round(rgb[0]*255)}, {round(rgb[1]*255)}, {round(rgb[2]*255)})"

