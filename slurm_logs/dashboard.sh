#!/bin/bash
# run this from the slurm_logs directory to get a simple tmux dashboard of the current run
# it shows: htop on the first node, a live view of the most recent slurm output, and a live view of the slurm queue

# Base name of the tmux session
BASE_SESSION_NAME="ARjl-dashboard"
SESSION_NAME=$BASE_SESSION_NAME

# Check if a session with the base name already exists
if tmux has-session -t $BASE_SESSION_NAME 2>/dev/null; then
    # Find a unique session name
    i=1
    while tmux has-session -t "${BASE_SESSION_NAME}_$i" 2>/dev/null; do
        ((i++))
    done
    SESSION_NAME="${BASE_SESSION_NAME}_$i"
fi

# Find the most recent SLURM output files
SLURM_OUT=$(ls -t ar*.out | head -n1)
SLURM_ERR=$(ls -t ar*.err | head -n1)

# Create a new tmux session
tmux new-session -d -s $SESSION_NAME

# Function to create a new pane and run a command
create_pane() {
    COMMAND=$1
    tmux split-window -t $SESSION_NAME
    tmux send-keys -t $SESSION_NAME "$COMMAND" C-m
}

tmux set -g pane-border-status top
tmux set -g pane-border-format "#{pane_index} #{pane_current_command}"

# tail the most recent SLURM output file
tmux send-keys -t $SESSION_NAME "tail -f \"\$(ls ar*.out | sort | tail -n1)\"" C-m
create_pane "ssh $(squeue -h -o '%N' -u $USER) -t htop -u $USER"
create_pane "while true; do clear; sq; sleep 10; done"

# Arrange the panes
tmux select-layout -t $SESSION_NAME main-vertical

# make the right-hand column a little wider
tmux resize-pane -t $SESSION_NAME -L 20

tmux attach-session -t $SESSION_NAME
