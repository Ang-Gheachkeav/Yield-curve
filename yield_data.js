// ==============================================================================
//  NBC should edit the following data set
// ------------------------------------------------------------------------------
//  target Fields: KHR and USD only, (others will be calculated automatically)
// ============================================================================== 
class Defaults {
	static dev_mode = 0; //0 for Production Mode
	static x_days = [7, 14, 28, 91, 182, 365, 730, 1095, 1825, 3650]; //Days NOT Years
	static asAt = "Nov 2025";
	static polyOrder = 2;  //2 for x^2, 3 for x^3
	static NSS_model = 2;  //1: Manual Data, 2: Auto Data, 3: Both
	static chart_data = {
		//Maintain the chart data here
		"KHR"	: {
			"name"	: "KHR",
			"color" : "rgba(255,140,64,.9)",
			"line"	: "solid",
			"data"	: [0.841710758377425, null, null, 0.876275787187839, 1.2902380952381, 3.4375, 4.27307692307692, 4.58891954022988, 5.25, 6.00],
		//	"data"	: [0.84, null, null, 0.88, 1.29, 3.44, 4.27, 4.59, 5.25, 6.00],
		},
		"InterBank in KHR"	:{
			"name"	: "Interbank in KHR",
			"type"	: "solid",
			"color"	: "rgba(156,217,107,.9)",
			"data"	: [2.87, 3.51, 3.76, 4.25, 4.58, 5.61, 5.85, 6.24, 6.86, 6.70]
		},
		"NSS_IBKHR"	:{
			"name"	: "NSS (Interbank KHR) Auto",
			"type"	: "solid",
			"color"	: "rgba(255, 255, 255, 0.9)",
			"data"	: [], //Auto Caluculation
		},
		"KHRNSS_MAN": {
			"name"	: "NSS (KHR) Manual",
			"color"	: "rgba(102,227,161,.9)",
			"line"	: "solid",
			"data"	: [0.43, 0.50, 0.63, 1.19, 1.88, 2.92, 4.12, 4.78, 5.31, 5.93],
		},
		"KHRNSS_AUTO": {
			"name"	: "NSS (KHR) Auto",
			"color"	: "rgba(236, 233, 76, 0.9)",
			"line"	: "dashed",
			"data"	: [], //Auto Caluculation
		},
		"USD": {
			"name"	: "USD",
			"color"	: "rgba(106,167,255,.9)",
			"line"	: "solid",
			"data"	: [null, null, null, 0.15, 0.63, null, null, null, null, null],
		},
		"USDNSS": {
			"name"	: "NSS (USD)",
			"color"	: "rgba(255,200,160,.9)",
			"line"	: "dashed",
			"data"	: [0.00, 0.02, 0.07, 0.28, 0.57, 1.12, null, null, null, null],
		},
		"KHRPOLY": {
			"name"	: "Polynomial (KHR)",
			"color"	: "rgba(180,100,220,.9)",
			"line"	: "dashed",
			"data"	: [], //Auto Caluculation
		},
	};
}
