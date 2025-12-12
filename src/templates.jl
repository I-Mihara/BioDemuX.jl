# Helper to generate SVG bar chart
function generate_histogram(data::Dict{T,Int}, title::String, x_label::String) where T
    if isempty(data)
        return "<p>No data for $title</p>"
    end

    sorted_keys = sort(collect(keys(data)))
    values = [data[k] for k in sorted_keys]
    max_val = maximum(values)



    width = 600
    height = 300
    margin_left = 60
    margin_bottom = 60
    bar_width = (width - margin_left) / length(values) * 0.8



    bar_color = "#3498db"

    svg = """<div class="chart-container"><h3>$title</h3><svg width="$width" height="$height" viewBox="0 0 $width $height">"""



    for i in 0:5
        y = height - margin_bottom - (i / 5) * (height - margin_bottom - 20)
        svg *= """<line x1="$margin_left" y1="$y" x2="$width" y2="$y" stroke="#eee" stroke-width="1" />"""
        # Y-axis labels
        label_val = round(Int, (i / 5) * max_val)
        svg *= """<text x="$(margin_left - 10)" y="$(y + 5)" font-size="10" text-anchor="end" fill="#7f8c8d">$label_val</text>"""
    end

    # Y-axis title
    svg *= """<text transform="rotate(-90)" x="$(-(height - margin_bottom)/2)" y="15" font-size="12" text-anchor="middle" fill="#2c3e50" font-weight="bold">Read Count</text>"""

    # Y-axis line
    svg *= """<line x1="$margin_left" y1="0" x2="$margin_left" y2="$(height - margin_bottom)" stroke="#ccc" />"""
    # X-axis line
    svg *= """<line x1="$margin_left" y1="$(height - margin_bottom)" x2="$width" y2="$(height - margin_bottom)" stroke="#ccc" />"""

    for (i, val) in enumerate(values)
        bar_height = (val / max_val) * (height - margin_bottom - 20)
        x = margin_left + (i - 1) * ((width - margin_left) / length(values)) + ((width - margin_left) / length(values) - bar_width) / 2
        y = height - margin_bottom - bar_height

        # Bar with hover effect (using simple title for now, CSS hover can be added)
        svg *= """<rect x="$x" y="$y" width="$bar_width" height="$bar_height" fill="$bar_color" class="bar"><title>$(sorted_keys[i]): $val</title></rect>"""

        # X-axis label
        if length(values) <= 20 || i % (div(length(values), 20) + 1) == 0
            label_x = x + bar_width / 2
            label_y = height - margin_bottom + 15
            svg *= """<text x="$label_x" y="$label_y" font-size="10" text-anchor="middle" fill="#7f8c8d">$(sorted_keys[i])</text>"""
        end
    end

    # X-axis title
    svg *= """<text x="$(width/2 + margin_left/2)" y="$(height - 10)" font-size="12" text-anchor="middle" fill="#2c3e50" font-weight="bold">$x_label</text>"""

    svg *= "</svg></div>"
    return svg
end

const HTML_REPORT_TEMPLATE_START = """
<!DOCTYPE html>
<html>
<head>
    <title>BioDemuX Report</title>
    <style>
        body { font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background-color: #f0f2f5; color: #1a1a1a; margin: 0; padding: 40px; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 4px 20px rgba(0,0,0,0.05); }
        h1 { color: #1a1a1a; border-bottom: 2px solid #f0f2f5; padding-bottom: 20px; margin-bottom: 30px; font-weight: 700; letter-spacing: -0.5px; }
        h2 { color: #2c3e50; margin-top: 40px; font-weight: 600; font-size: 1.5em; }
        h3 { color: #57606a; font-size: 1em; font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; }
        
        .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 24px; margin-bottom: 40px; }
        .summary-card { background: white; padding: 24px; border-radius: 12px; border: 1px solid #e1e4e8; box-shadow: 0 2px 8px rgba(0,0,0,0.02); transition: transform 0.2s; }
        .summary-card:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.05); }
        .summary-card h3 { margin: 0 0 12px 0; font-size: 0.85em; color: #6e7781; }
        .summary-card .value { font-size: 2.5em; font-weight: 700; color: #1a1a1a; line-height: 1; margin-bottom: 8px; }
        .summary-card .percentage { font-size: 1em; color: #57606a; font-weight: 500; }
        
        table { border-collapse: separate; border-spacing: 0; width: 100%; margin-top: 20px; border-radius: 8px; overflow: hidden; border: 1px solid #e1e4e8; }
        th, td { padding: 16px 20px; text-align: left; border-bottom: 1px solid #e1e4e8; }
        th { background-color: #f6f8fa; font-weight: 600; color: #24292f; font-size: 0.9em; text-transform: uppercase; letter-spacing: 0.5px; }
        tr:last-child td { border-bottom: none; }
        tr:hover td { background-color: #f8f9fa; }
        
        .progress-bar-container { width: 100%; background-color: #e1e4e8; border-radius: 100px; height: 8px; position: relative; overflow: hidden; }
        .progress-bar { height: 100%; background-color: #3498db; border-radius: 100px; }
        .progress-label { display: none; } /* Hide label inside bar, show count in table cell */
        
        .chart-container { margin-top: 24px; padding: 24px; background: white; border: 1px solid #e1e4e8; border-radius: 12px; overflow-x: auto; }
        .bar:hover { opacity: 0.8; cursor: pointer; }
        
        .run-info { background: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 30px; border: 1px solid #e1e4e8; }
        .run-info h3 { margin-top: 0; color: #2c3e50; }
        .run-info ul { list-style: none; padding: 0; margin: 0; }
        .run-info li { margin-bottom: 8px; font-size: 0.95em; color: #57606a; }
        .run-info li strong { color: #24292f; font-weight: 600; margin-right: 8px; }

        .tabs { display: flex; border-bottom: 1px solid #e1e4e8; margin-top: 40px; gap: 8px; }
        .tab { padding: 12px 24px; cursor: pointer; background: transparent; border: none; border-bottom: 2px solid transparent; font-weight: 500; color: #57606a; transition: all 0.2s; }
        .tab:hover { color: #24292f; }
        .tab.active { border-bottom: 2px solid #3498db; color: #3498db; font-weight: 600; }
        .tab-content { display: none; padding: 30px 0; animation: fadeIn 0.3s ease; }
        .tab-content.active { display: block; }
        
        @keyframes fadeIn { from { opacity: 0; transform: translateY(10px); } to { opacity: 1; transform: translateY(0); } }

        select { padding: 8px 12px; border: 1px solid #e1e4e8; border-radius: 6px; font-size: 14px; color: #24292f; background-color: white; min-width: 200px; }
    </style>
    <script>
        function openTab(evt, tabName) {
            var tabButton = evt.currentTarget;
            var tabContainer = tabButton.parentElement;
            var allTabs = tabContainer.getElementsByClassName("tab");
            
            for (var i = 0; i < allTabs.length; i++) {
                allTabs[i].classList.remove("active");
            }
            
            tabButton.classList.add("active");
            
            var parts = tabName.split("_");
            var id = parts[parts.length-1];
            
            document.getElementById("Global_" + id).style.display = "none";
            document.getElementById("PerBarcode_" + id).style.display = "none";
            document.getElementById("Global_" + id).classList.remove("active");
            document.getElementById("PerBarcode_" + id).classList.remove("active");
            
            document.getElementById(tabName).style.display = "block";
            document.getElementById(tabName).classList.add("active");
        }

        function showBarcodeStats(select, runId) {
            var bcId = select.value;
            var statsDivs = document.getElementsByClassName("barcode-stats-" + runId);
            for (var i = 0; i < statsDivs.length; i++) {
                statsDivs[i].style.display = "none";
            }
            if (bcId) {
                document.getElementById("stats-" + bcId + "-" + runId).style.display = "block";
            }
        }
    </script>
</head>
<body>
    <div class="container">
        <h1>BioDemuX Report</h1>
"""

const HTML_REPORT_TEMPLATE_END = """
    </div>
</body>
</html>
"""
