"""Terminal presentation for the MetaQuest command-line interface."""

import sys
import time
import subprocess
import threading
import re
from pathlib import Path
from datetime import datetime
from contextlib import contextmanager
from typing import Optional, Dict, List, Any, Tuple

class Colors:
    """ANSI color codes for terminal output"""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'
    
    @classmethod
    def disable(cls):
        """Disable colors for non-terminal outputs"""
        cls.HEADER = cls.BLUE = cls.CYAN = cls.GREEN = ''
        cls.YELLOW = cls.RED = cls.BOLD = cls.DIM = cls.UNDERLINE = cls.END = ''


class Spinner:
    """Animated spinner for long-running operations - Thread-safe version"""
    def __init__(self, message: str = "Processing"):
        self.message = message
        self._spinning = False  # Changed to private
        self._lock = threading.Lock()  # NEW: Thread safety
        self.thread = None
        self.frames = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
        self._frame_idx = 0  # Changed to private
        self._stream = sys.stderr
    
    @property
    def spinning(self):
        """Thread-safe getter for spinning state"""
        with self._lock:
            return self._spinning
    
    @spinning.setter
    def spinning(self, value):
        """Thread-safe setter for spinning state"""
        with self._lock:
            self._spinning = value
    
    def _spin(self):
        """Spinner animation loop with thread-safe access"""
        while self.spinning:  # Now uses thread-safe property
            try:
                if self._stream and not self._stream.closed:
                    # Protect frame_idx from concurrent access
                    with self._lock:
                        frame = self.frames[self._frame_idx % len(self.frames)]
                        self._frame_idx += 1
                    
                    self._stream.write(f'\r  {Colors.CYAN}[{frame}]{Colors.END} {self.message}...')
                    self._stream.flush()
                else:
                    break
            except (ValueError, AttributeError, OSError):
                break
            time.sleep(0.1)
    
    def start(self):
        """Start the spinner"""
        if not self._stream or self._stream.closed:
            return
            
        self.spinning = True  # Uses thread-safe setter
        self.thread = threading.Thread(target=self._spin, daemon=True)
        try:
            self.thread.start()
        except RuntimeError:
            self.spinning = False
    
    def stop(self, final_message: Optional[str] = None):
        """Stop the spinner gracefully with thread-safe cleanup"""
        self.spinning = False  # Uses thread-safe setter
        
        # Wait for thread to finish
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout=0.5)
        
        # Clear the spinner line with exclusive access
        try:
            if self._stream and not self._stream.closed:
                with self._lock:  # Ensure no concurrent writes
                    if final_message:
                        self._stream.write(f'\r  {final_message}\n')
                    else:
                        self._stream.write('\r' + ' ' * (len(self.message) + 20) + '\r')
                    self._stream.flush()
        except (ValueError, AttributeError, OSError):
            pass

class TableFormatter:
    """Beautiful table formatting with box-drawing characters"""
    
    @staticmethod
    def format_table(title: str, sections: List[Dict[str, Any]], width: int = 73) -> str:
        """Format a beautiful table with sections"""
        lines = []
        
        # Top border with title
        lines.append(f"┌{'─' * (width - 2)}┐")
        title_padding = width - len(title) - 3
        lines.append(f"│ {Colors.BOLD}{title}{Colors.END}{' ' * title_padding}│")
        
        for i, section in enumerate(sections):
            # Section divider
            lines.append(f"├{'─' * (width - 2)}┤")
            
            # Section header (if provided)
            if 'header' in section and section['header']:
                header = section['header']
                header_padding = width - len(header) - 3
                lines.append(f"│ {Colors.CYAN}{header}{Colors.END}{' ' * header_padding}│")
                lines.append(f"├{'─' * (width - 2)}┤")
            
            # Section rows
            for key, value in section.get('rows', {}).items():
                # Skip rows with None or 0/100 values unless they're meaningful
                if value is None or (isinstance(value, str) and value == "0/100"):
                    continue
                    
                key_display = f"  {key}"
                value_display = str(value)
                
                # Calculate spacing
                key_visible = len(key_display)
                value_visible = len(value_display)
                space_needed = width - key_visible - value_visible - 3
                
                if space_needed < 1:
                    space_needed = 1
                
                lines.append(f"│{key_display}{' ' * space_needed}{value_display} │")
        
        # Bottom border
        lines.append(f"└{'─' * (width - 2)}┘")
        
        return '\n'.join(lines)
    
    @staticmethod
    def format_visual_bar(percentage: float, width: int = 10, label: str = "") -> str:
        """Create a visual percentage bar"""
        filled = int(width * (percentage / 100))
        bar = '█' * filled + '░' * (width - filled)
        
        # Color based on percentage
        if percentage >= 90:
            color = Colors.GREEN
        elif percentage >= 70:
            color = Colors.YELLOW
        else:
            color = Colors.RED
        
        result = f"{bar} {color}{percentage:5.1f}%{Colors.END}"
        if label:
            result += f" {label}"
        
        return result


class OutputFormatter:
    """Enhanced output formatter with better verbosity control"""
    
    SILENT = 0      # No output except errors
    STANDARD = 1    # Normal user output
    VERBOSE = 2     # Additional details
    DEBUG = 3       # Everything including debug

    def __init__(
        self,
        verbosity: str = 'standard',
        log_file: Optional[Path] = None,
        *,
        append_log: bool = False,
    ):
        self.verbosity_map = {
            'silent': self.SILENT,
            'standard': self.STANDARD,
            'verbose': self.VERBOSE,
            'debug': self.DEBUG
        }
        self.verbosity = self.verbosity_map.get(verbosity, self.STANDARD)
        self.colors_enabled = sys.stderr.isatty()
        self.log_file = log_file
        self.start_time = time.time()
        self._active_spinner = None  # Track current spinner for debug output pausing
        
        if not self.colors_enabled:
            Colors.disable()
        
        if self.log_file:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            self.log_handle = open(
                self.log_file,
                'a' if append_log else 'w',
                encoding='utf-8',
            )
            if append_log and self.log_file.stat().st_size:
                self.log_handle.write("\n=== MetaQuest resumed run ===\n")
                self.log_handle.flush()
        else:
            self.log_handle = None

    def __del__(self):
        if hasattr(self, 'log_handle') and self.log_handle:
            self.log_handle.close()

    def _log(self, message: str):
        """Write to log file if enabled"""
        if self.log_handle:
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            clean_msg = re.sub(r'\x1b\[[0-9;]*m', '', message)
            self.log_handle.write(f"[{timestamp}] {clean_msg}\n")
            self.log_handle.flush()

    def _print(self, message: str, min_verbosity: int = STANDARD):
        """Internal print with verbosity check"""
        if self.verbosity >= min_verbosity:
            print(message, file=sys.stderr)
        self._log(message)

    # ========================================================================
    # SECTION HEADERS
    # ========================================================================
    
    def banner(self, title: str, version: str, tagline: str):
        """Print a compact application header."""
        if self.verbosity >= self.STANDARD:
            self._print(
                f"{Colors.BOLD}{Colors.CYAN}{title}{Colors.END} "
                f"{Colors.DIM}{version} · {tagline}{Colors.END}"
            )
            self._log(f"{title} v{version} - {tagline}")

    def section_header(self, title: str):
        """Print a restrained subsection label."""
        if self.verbosity >= self.STANDARD:
            self._log("")
            self._log(title)
            self._print(f"\n{Colors.BOLD}{title.title()}{Colors.END}")

    # ========================================================================
    # STATUS MESSAGES - Simplified tree structure
    # ========================================================================
    
    def substep(self, message: str):
        """Print sub-step with deeper tree level"""
        if self.verbosity >= self.VERBOSE:
            self._print(f"        • {message}")

    def success(self, message: str, style: Optional[str] = None):
        """Print success message"""
        if self.verbosity >= self.STANDARD:
            if style == 'bold':
                message = f"{Colors.BOLD}{message}{Colors.END}"
            self._print(f"  {Colors.GREEN}✓{Colors.END} {message}")

    def info(self, message: str, indent: int = 1):
        """Print info message"""
        if self.verbosity >= self.STANDARD:
            prefix = "  " if indent == 1 else "    "
            self._print(f"{prefix}{Colors.DIM}·{Colors.END} {message}")

    def warning(self, message: str):
        """Print warning message"""
        if self.verbosity >= self.STANDARD:
            self._print(f"  {Colors.YELLOW}!{Colors.END} {message}")

    def error(self, message: str, solutions: Optional[List[str]] = None, doc_link: Optional[str] = None):
        """Print error message with solutions"""
        if self.verbosity >= self.SILENT:  # Always show errors
            self._print(
                f"\n     {Colors.RED}✗ {message}{Colors.END}",
                min_verbosity=self.SILENT,
            )
            if solutions:
                self._print(
                    f"\n     {Colors.BOLD}💡 Solutions:{Colors.END}",
                    min_verbosity=self.SILENT,
                )
                for solution in solutions:
                    self._print(
                        f"        • {solution}",
                        min_verbosity=self.SILENT,
                    )
            if doc_link:
                self._print(
                    f"\n     {Colors.CYAN}📚 Documentation:{Colors.END} {doc_link}",
                    min_verbosity=self.SILENT,
                )
            print(file=sys.stderr)

    # ========================================================================
    # PROGRESS TRACKING
    # ========================================================================
    
    @contextmanager
    def spinner(self, message: str):
        """Context manager for spinner animation"""
        if self.verbosity >= self.STANDARD and self.colors_enabled:
            spinner = Spinner(message)
            self._active_spinner = spinner
            spinner.start()
            try:
                yield
            finally:
                spinner.stop()
                self._active_spinner = None
        else:
            # In silent mode, just log without printing
            if self.verbosity >= self.VERBOSE:
                self.info(message)
            yield

    # ========================================================================
    # SUBPROCESS EXECUTION with better output control
    # ========================================================================
    
    def run_subprocess(self, cmd: List[str], operation_name: str, 
                      capture_output: bool = True, show_command: bool = False) -> Tuple[int, str, str]:
        """Run subprocess with controlled output. All tool output is structured and timestamped."""
        if show_command or self.verbosity >= self.DEBUG:
            self.debug(f"Command: {' '.join(cmd)}")
        
        if capture_output:
            if self.verbosity >= self.DEBUG:
                # DEBUG mode: capture output, print through structured block (no raw sys.stderr writes)
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,
                    encoding='utf-8',
                    errors='replace'
                )
                
                stdout_lines = []
                stderr_lines = []
                
                def read_debug_stream(stream, line_list, tag):
                    """Collect lines without writing to sys.stderr/stdout directly."""
                    for line in stream:
                        line_list.append(line)
                
                t1 = threading.Thread(target=read_debug_stream, args=(process.stdout, stdout_lines, "STDOUT"), daemon=True)
                t2 = threading.Thread(target=read_debug_stream, args=(process.stderr, stderr_lines, "STDERR"), daemon=True)
                t1.start()
                t2.start()
                
                return_code = process.wait()
                t1.join(timeout=2)
                t2.join(timeout=2)
                
                stdout_str = "".join(stdout_lines)
                stderr_str = "".join(stderr_lines)
                
                # Print structured block for STDOUT, then STDERR
                if stdout_str.strip():
                    self._debug_block_header(operation_name, "STDOUT")
                    for line in stdout_str.splitlines():
                        self._debug_tool_line(line, f"STDOUT:{operation_name}")
                    self._debug_block_footer(operation_name, "STDOUT")
                else:
                    # Still log even if nothing to print
                    if stdout_str:
                        self._log(f"[STDOUT:{operation_name}] --- begin ---")
                        for line in stdout_str.splitlines():
                            self._log(f"[STDOUT:{operation_name}] {line}")
                        self._log(f"[STDOUT:{operation_name}] --- end ---")

                if stderr_str.strip():
                    self._debug_block_header(operation_name, "STDERR")
                    for line in stderr_str.splitlines():
                        self._debug_tool_line(line, f"STDERR:{operation_name}")
                    self._debug_block_footer(operation_name, "STDERR")
                else:
                    if stderr_str:
                        self._log(f"[STDERR:{operation_name}] --- begin ---")
                        for line in stderr_str.splitlines():
                            self._log(f"[STDERR:{operation_name}] {line}")
                        self._log(f"[STDERR:{operation_name}] --- end ---")
                
                return return_code, stdout_str, stderr_str
            else:
                # Non-debug: fully suppress terminal output, log everything
                result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True, encoding='utf-8')
                if result.stdout:
                    self._log(f"[STDOUT:{operation_name}] --- begin ---")
                    for line in result.stdout.splitlines():
                        self._log(f"[STDOUT:{operation_name}] {line}")
                    self._log(f"[STDOUT:{operation_name}] --- end ---")
                if result.stderr:
                    self._log(f"[STDERR:{operation_name}] --- begin ---")
                    for line in result.stderr.splitlines():
                        self._log(f"[STDERR:{operation_name}] {line}")
                    self._log(f"[STDERR:{operation_name}] --- end ---")
                return result.returncode, result.stdout, result.stderr
        else:
            # capture_output=False: always capture so we can log/display cleanly
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, encoding='utf-8', errors='replace')
            stdout_str = result.stdout or ""
            stderr_str = result.stderr or ""
            if stdout_str.strip():
                self._debug_block_header(operation_name, "STDOUT")
                for line in stdout_str.splitlines():
                    self._debug_tool_line(line, f"STDOUT:{operation_name}")
                self._debug_block_footer(operation_name, "STDOUT")
            elif stdout_str:
                self._log(f"[STDOUT:{operation_name}] --- begin ---")
                for line in stdout_str.splitlines():
                    self._log(f"[STDOUT:{operation_name}] {line}")
                self._log(f"[STDOUT:{operation_name}] --- end ---")
            if stderr_str.strip():
                self._debug_block_header(operation_name, "STDERR")
                for line in stderr_str.splitlines():
                    self._debug_tool_line(line, f"STDERR:{operation_name}")
                self._debug_block_footer(operation_name, "STDERR")
            elif stderr_str:
                self._log(f"[STDERR:{operation_name}] --- begin ---")
                for line in stderr_str.splitlines():
                    self._log(f"[STDERR:{operation_name}] {line}")
                self._log(f"[STDERR:{operation_name}] --- end ---")
            return result.returncode, stdout_str, stderr_str
            
    # ========================================================================
    # UTILITY FUNCTIONS
    # ========================================================================
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration"""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.0f}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours}h {minutes}m"

    def debug(self, message: str):
        """Print debug message"""
        if self.verbosity >= self.DEBUG:
            self._print(f"     {Colors.DIM}[DEBUG] {message}{Colors.END}")

    def _debug_tool_line(self, line: str, tag: str):
        """
        Print a single tool output line in debug mode.
        Goes to terminal as a structured dim line — never writes raw to sys.stderr/stdout.
        Also logs with timestamp via _log().
        """
        stripped = line.rstrip()
        # Print to terminal only — clean, indented, dim
        if self.verbosity >= self.DEBUG:
            print(f"  {Colors.DIM}│ {stripped}{Colors.END}", file=sys.stderr)
        # Always log with timestamp (separate from terminal print)
        self._log(f"[{tag}] {stripped}")

    def _debug_block_header(self, operation_name: str, tag: str):
        """Print a clean visual block header for tool output in debug mode."""
        if self.verbosity >= self.DEBUG:
            # Stop + clear the spinner line before printing (prevents interleaving)
            if self._active_spinner and self._active_spinner.spinning:
                self._active_spinner.spinning = False
                if self._active_spinner.thread:
                    self._active_spinner.thread.join(timeout=0.3)
                # Clear the spinner line
                sys.stderr.write('\r\033[K')
                sys.stderr.flush()
            label = f" {tag}: {operation_name} "
            bar = "─" * max(0, (73 - len(label)) // 2)
            print(f"\n  {Colors.DIM}{bar}{label}{bar}{Colors.END}", file=sys.stderr)
        self._log(f"[{tag}:{operation_name}] --- begin ---")

    def _debug_block_footer(self, operation_name: str, tag: str):
        """Print a clean visual block footer for tool output in debug mode."""
        if self.verbosity >= self.DEBUG:
            print(f"  {Colors.DIM}{'─' * 73}{Colors.END}\n", file=sys.stderr)
        self._log(f"[{tag}:{operation_name}] --- end ---")

    def result(self, metrics: Dict[str, Any], indent: int = 1):
        """Display aligned key/value metrics."""
        if self.verbosity >= self.STANDARD:
            prefix = "    " if indent == 2 else "  "
            width = max((len(str(key)) for key in metrics), default=0)
            for key, value in metrics.items():
                self._print(
                    f"{prefix}{Colors.DIM}{str(key).ljust(width)}{Colors.END}  {value}"
                )


# ========================================================================
# GLOBAL FORMATTER INSTANCE
# ========================================================================
_formatter_instance = None

def get_formatter() -> OutputFormatter:
    global _formatter_instance
    if _formatter_instance is None:
        _formatter_instance = OutputFormatter()
    return _formatter_instance

def set_formatter(formatter: OutputFormatter):
    global _formatter_instance
    _formatter_instance = formatter
