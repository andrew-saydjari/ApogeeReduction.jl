# This script generates a dashboard of diagnostic plots for a given MJD.
# Example usage (from the root of the ApogeeReduction.jl repository):
# julia +1.11.0 --project scripts/run/generate_dashboard.jl --mjd 60835 --outdir ../outdir/

using ArgParse, Glob, SlackThreads
using ApogeeReduction # for the slack environment variables

########################################################
# Categorize the plots
# edit here to categorize new plot types
########################################################
# Define the main categories and their patterns in the desired order
# Each of these is a section name and a regex pattern for what should be in that section
CATEGORY_PATTERNS = [
    ("ar1Dunical", r"ar1Dunical_.*\.png$"),
    ("ar1Duni", r"ar1Duni_.*\.png$"),
    ("wave", r".*wave_.*\.png$"),
    ("night wave", r"night.*wave_.*\.png$"),
    ("skyPeakResiduals", r"skyPeakResiduals_.*\.png$"),
    ("ar1D", r"ar1D_.*\.png$"),
    ("ar2Dresidualscal", r"ar2Dresidualscal_.*\.png$"),
    ("fpi", r"fpi.*\.png$"),
    ("Parameters", r".*Params_.*\.png$")
]

# Define the order in which object types should appear
# Leave empty to use alphabetical order
# the empty string is for files that don't have an object type
const OBJECT_TYPES = ["object", "domeflat", "internalflat", "quartzflat", "dark", "arclamp", ""]
const OBJECT_TYPE_ORDER = Dict(OBJECT_TYPES .=> 1:length(OBJECT_TYPES))

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

function extract_object_type(filename)
    # Extract object type from filename like "something_something_something_OBJECTTYPE.png"
    idx = findlast('_', basename(filename))
    obj_type = basename(filename)[(idx + 1):(end - 4)] # cut off the .png extension
    if all(islowercase(c) for c in obj_type)
        obj_type
    else
        ""
    end
end

function get_plot_categories(plot_files)
    categories = []
    for (category, pattern) in CATEGORY_PATTERNS
        matching_files = filter(f -> occursin(pattern, f), plot_files)
        isempty(matching_files) && continue

        # Group by object type
        object_groups = Dict{String, Vector{String}}()
        for file in matching_files
            obj_type = extract_object_type(file)
            haskey(object_groups, obj_type) || (object_groups[obj_type] = String[])
            push!(object_groups[obj_type], file)
        end

        # Sort object types according to OBJECT_TYPE_ORDER
        obj_types = sort(collect(keys(object_groups)),
            by = obj_type -> get(OBJECT_TYPE_ORDER, obj_type, 1000))

        # Create ordered list of (object_type, files) tuples
        object_sections = [(obj_type, sort(object_groups[obj_type])) for obj_type in obj_types]
        push!(categories, (category, object_sections))
    end

    # put files that don't match any category in "other" and print a warning
    other_files = setdiff(plot_files,
        [f for (_, obj_sections) in categories for (_, files) in obj_sections for f in files])
    if !isempty(other_files)
        println("Warning: $(length(other_files)) files don't match any category. (Adding to 'other' category)")
        # Group other files by object type too
        push!(categories, ("other", [("", other_files)]))
    end

    categories
end

function generate_html(mjd, outdir, categories)
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>ApogeeReduction.jl diagnostic plots - MJD $mjd</title>
        <style>
            :root {
                --primary-color: #00ff00;
                --background-color: #000000;
                --secondary-bg: #1a1a1a;
                --text-color: #ffffff;
                --secondary-text: #cccccc;
                --border-color: #333;
            }

            body {
                background-color: var(--background-color);
                color: var(--text-color);
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
                color: var(--primary-color);
                border-bottom: 1px solid var(--primary-color);
                padding-bottom: 10px;
            }
            .category h3 {
                color: var(--secondary-text);
                margin-top: 20px;
                margin-bottom: 10px;
                font-size: 1.1em;
            }
            .plot-grid {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
                gap: 20px;
                margin-top: 20px;
            }
            .plot-item {
                background-color: var(--secondary-bg);
                padding: 10px;
                border-radius: 5px;
                cursor: pointer;
            }
            .plot-item img {
                width: 100%;
                height: auto;
                border-radius: 3px;
                transition: transform 0.2s;
                display: block;
                object-fit: contain;
            }
            .plot-title {
                margin-top: 10px;
                font-size: 0.9em;
                color: var(--secondary-text);
                word-break: break-word;
            }
            .modal {
                display: none;
                position: fixed;
                z-index: 1000;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background-color: rgba(0, 0, 0, 0.7);
                padding: 20px;
                box-sizing: border-box;
            }
            .modal-content {
                max-width: 90%;
                max-height: 80vh;
                margin: auto;
                display: block;
            }
            .nav-button {
                position: fixed;
                bottom: 20px;
                transform: none;
                background-color: rgba(0, 0, 0, 0.7);
                color: var(--text-color);
                border: 2px solid var(--primary-color);
                border-radius: 50%;
                width: 80px;
                height: 80px;
                font-size: 24px;
                cursor: pointer;
                display: flex;
                align-items: center;
                justify-content: center;
                z-index: 1001;
            }
            .nav-button:hover {
                background-color: var(--primary-color);
                color: var(--background-color);
            }
            .nav-button.prev {
                left: 20px;
            }
            .nav-button.next {
                right: 20px;
            }
            .modal-filename {
                color: var(--text-color);
                text-align: center;
                margin-top: 15px;
                font-size: 1.2em;
            }
            .modal-hint {
                color: #888888;
                text-align: center;
                margin-top: 10px;
                font-size: 0.9em;
                font-style: italic;
            }
            .close {
                position: absolute;
                top: 15px;
                right: 35px;
                color: #f1f1f1;
                font-size: 40px;
                font-weight: bold;
                cursor: pointer;
            }
            .nav {
                position: sticky;
                top: 0;
                background-color: var(--background-color);
                padding: 10px 0;
                margin-bottom: 20px;
                border-bottom: 1px solid var(--border-color);
            }
            .nav a {
                color: var(--primary-color);
                text-decoration: none;
                margin-right: 20px;
                padding: 5px 10px;
                border-radius: 3px;
            }
            .nav a:hover {
                background-color: var(--secondary-bg);
            }
            .nav a.active {
                background-color: var(--secondary-bg);
                border: 1px solid var(--primary-color);
            }
            .log-links a {
                color: var(--primary-color);
                text-decoration: none;
                margin-right: 20px;
                padding: 5px 10px;
            }
        </style>
    </head>
    <body>
        <div class="container">
                <h1>ApogeeReduction.jl diagnostic plots - MJD $mjd</h1>
                <div class="log-links">
                    <a href="https://data.sdss5.org/sas/sdsswork/mwm/sandbox/airflow-ApogeeReduction.jl/daily/ApogeeReduction.jl/metadata/observing_log_viewer/?sjd=$mjd&site=apo" target="_blank">APO Observing Log</a>
                    <a href="https://data.sdss5.org/sas/sdsswork/mwm/sandbox/airflow-ApogeeReduction.jl/daily/ApogeeReduction.jl/metadata/observing_log_viewer/?sjd=$mjd&site=lco" target="_blank">LCO Observing Log</a>
                </div>
            <div class="nav">
    """

    # Add navigation links
    for (category, _) in categories
        category_id = replace(lowercase(category), " " => "-")
        html *= """<a href="#$category_id">$category</a>\n"""
    end
    html *= "</div>"

    # Add content sections
    for (category, object_sections) in categories
        category_id = replace(lowercase(category), " " => "-")
        html *= """
            <div class="category" id="$category_id">
                <h2>$category</h2>
        """

        for (object_type, files) in object_sections
            html *= """
                <h3>$object_type</h3>
                <div class="plot-grid">
            """
            for file in files
                filename = basename(file)
                html *= """
                        <div class="plot-item" onclick="openModal(this)">
                            <img src="$filename" alt="$filename">
                            <div class="plot-title">$filename</div>
                        </div>
                """
            end
            html *= "</div>"
        end
        html *= "</div>"
    end

    html *= """
        </div>
        <div id="imageModal" class="modal" onclick="closeModal()">
            <span class="close">&times;</span>
            <img class="modal-content" id="modalImage">
            <button class="nav-button prev" onclick="navigateImage(-1, event)">&larr;</button>
            <button class="nav-button next" onclick="navigateImage(1, event)">&rarr;</button>
            <div class="modal-filename" id="modalFilename"></div>
            <div class="modal-hint">Use &larr; and &rarr; arrow keys or click the buttons to navigate between plots</div>
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

            // Modal functionality
            let currentIndex = 0;
            let allPlotItems = [];

            function openModal(element) {
                const modal = document.getElementById('imageModal');
                const modalImg = document.getElementById('modalImage');
                const modalFilename = document.getElementById('modalFilename');

                // Get all plot items if we haven't already
                if (allPlotItems.length === 0) {
                    allPlotItems = Array.from(document.querySelectorAll('.plot-item'));
                }

                // Find the current index
                currentIndex = allPlotItems.indexOf(element);

                modal.style.display = "block";
                modalImg.src = element.querySelector('img').src;
                modalFilename.textContent = element.querySelector('.plot-title').textContent;
            }

            function closeModal() {
                document.getElementById('imageModal').style.display = "none";
            }

            function navigateImage(direction, event) {
                event.stopPropagation();
                currentIndex = (currentIndex + direction + allPlotItems.length) % allPlotItems.length;
                openModal(allPlotItems[currentIndex]);
            }

            // Close modal when clicking the close button
            document.querySelector('.close').onclick = function(event) {
                event.stopPropagation();
                closeModal();
            }

            // Add keyboard navigation
            document.addEventListener('keydown', function(event) {
                const modal = document.getElementById('imageModal');
                if (modal.style.display === "block") {
                    if (event.key === "ArrowLeft" || event.key === "Left") {
                        event.preventDefault();
                        currentIndex = (currentIndex - 1 + allPlotItems.length) % allPlotItems.length;
                        openModal(allPlotItems[currentIndex]);
                    } else if (event.key === "ArrowRight" || event.key === "Right") {
                        event.preventDefault();
                        currentIndex = (currentIndex + 1) % allPlotItems.length;
                        openModal(allPlotItems[currentIndex]);
                    }
                }
            });
        </script>
    </body>
    </html>
    """

    return html
end

function main()
    parg = parse_commandline()
    mjd = parg["mjd"]
    outdir = parg["outdir"]

    plots_dir = joinpath(outdir, "plots", string(mjd))
    plot_files = glob("*.png", plots_dir)
    categories = get_plot_categories(plot_files)
    html = generate_html(mjd, outdir, categories)
    output_file = joinpath(plots_dir, "dashboard.html")
    write(output_file, html)

    println("Dashboard generated at: $output_file")

    # post a link to slack
    # the link to the dashboard is relative to the current directory
    dir = abspath(joinpath(pwd(), parg["outdir"], "plots", string(mjd), "dashboard.html"))
    url = "<"*replace(replace(replace(dir, 
        "/uufs/chpc.utah.edu/common/home/sdss42/" => "https://data.sdss5.org/sas/"),
        "/mnt/ceph/users/sdssv/work/asaydjari/" => "https://users.flatironinstitute.org/~asaydjari/$(ENV["SLACK_TOKEN"])/sdsswork/"),
        "/mnt/ceph/users/sdssv/work/" => "https://users.flatironinstitute.org/~asaydjari/$(ENV["SLACK_TOKEN"])/"
    )*"|here>"
    thread = SlackThread()
    msg = "The reduction plots dashboard for SJD $(mjd) has been generated and is available at: $url"
    thread(msg)
    println(msg)
end

main()
