using ArgParse, Glob

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--mjd"
        required = true
        help = "MJD to generate dashboard for"
        arg_type = Int
        "--outdir"
        help = "Directory where data is stored"
        arg_type = String
        default = "../outdir/"
    end
    return parse_args(s)
end

function get_plot_categories(plot_files)
    categories = Dict{String, Vector{String}}()

    # Define the main categories and their patterns
    category_patterns = Dict(
        "Final Calibrated Spectra" => r"ar1Dunical_.*\.png$",
        "Uncalibrated Spectra" => r"ar1Duni_.*\.png$",
        "2D Residuals" => r"ar2Dresidualscal_.*\.png$",
        "Wavelength Solution" => r".*wave_.*\.png$",
        "Parameters" => r".*Params_.*\.png$"
    )

    # Sort files into categories
    for (category, pattern) in category_patterns
        matching_files = filter(f -> occursin(pattern, f), plot_files)
        if !isempty(matching_files)
            categories[category] = sort(matching_files)
        end
    end

    return categories
end

function generate_html(mjd, outdir, categories)
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>APOGEE Reduction Dashboard - MJD $mjd</title>
        <style>
            body {
                background-color: #000000;
                color: #ffffff;
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 20px;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
            }
            .category {
                margin-bottom: 30px;
            }
            .category h2 {
                color: #00ff00;
                border-bottom: 1px solid #00ff00;
                padding-bottom: 10px;
            }
            .plot-grid {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
                gap: 20px;
                margin-top: 20px;
            }
            .plot-item {
                background-color: #1a1a1a;
                padding: 10px;
                border-radius: 5px;
            }
            .plot-item img {
                width: 100%;
                height: auto;
                border-radius: 3px;
            }
            .plot-title {
                margin-top: 10px;
                font-size: 0.9em;
                color: #cccccc;
                word-break: break-word;
            }
            .nav {
                position: sticky;
                top: 0;
                background-color: #000000;
                padding: 10px 0;
                margin-bottom: 20px;
                border-bottom: 1px solid #333;
            }
            .nav a {
                color: #00ff00;
                text-decoration: none;
                margin-right: 20px;
                padding: 5px 10px;
                border-radius: 3px;
            }
            .nav a:hover {
                background-color: #1a1a1a;
            }
            .nav a.active {
                background-color: #1a1a1a;
                border: 1px solid #00ff00;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>APOGEE Reduction Dashboard - MJD $mjd</h1>

            <div class="nav">
    """

    # Add navigation links
    for (category, _) in categories
        category_id = replace(lowercase(category), " " => "-")
        html *= """<a href="#$category_id">$category</a>\n"""
    end

    html *= """
            </div>
    """

    # Add content sections
    for (category, files) in categories
        category_id = replace(lowercase(category), " " => "-")
        html *= """
            <div class="category" id="$category_id">
                <h2>$category</h2>
                <div class="plot-grid">
        """

        for file in files
            filename = basename(file)
            html *= """
                    <div class="plot-item">
                        <img src="$filename" alt="$filename">
                        <div class="plot-title">$filename</div>
                    </div>
            """
        end

        html *= """
                </div>
            </div>
        """
    end

    html *= """
        </div>
        <script>
            // Add active class to current section in nav
            const sections = document.querySelectorAll('.category');
            const navLinks = document.querySelectorAll('.nav a');

            window.addEventListener('scroll', () => {
                let current = '';
                sections.forEach(section => {
                    const sectionTop = section.offsetTop;
                    if (pageYOffset >= sectionTop - 60) {
                        current = section.getAttribute('id');
                    }
                });

                navLinks.forEach(link => {
                    link.classList.remove('active');
                    if (link.getAttribute('href').substring(1) === current) {
                        link.classList.add('active');
                    }
                });
            });
        </script>
    </body>
    </html>
    """

    return html
end

function main()
    args = parse_commandline()
    mjd = args["mjd"]
    outdir = args["outdir"]

    # Get the plots directory for this MJD
    plots_dir = joinpath(outdir, "plots", string(mjd))

    # Get all PNG files in the directory
    plot_files = glob("*.png", plots_dir)

    # Categorize the plots
    categories = get_plot_categories(plot_files)

    # Generate the HTML
    html = generate_html(mjd, outdir, categories)

    # Write the HTML file
    output_file = joinpath(plots_dir, "dashboard.html")
    write(output_file, html)

    println("Dashboard generated at: $output_file")
end

main()
